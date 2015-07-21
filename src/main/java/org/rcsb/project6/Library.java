package org.rcsb.project6;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.hadoop.io.Demo;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralAlignment.SuperPositionQCP;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;

import scala.Tuple2;
import scala.Tuple3;

/**
 * 
 * This class creates a library of unique fragments categorized by length.
 * Output PDBID as a string, length as a string, and the points as a point3d[].
 * 
 * @author Grant Summers
 *
 */

public class Library
{
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3;	
	private static int length = 8;
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException
	{
		 // arguments
		 String path = "/Users/grantsummers/Desktop/School/Internship/main/QuickStructureSearch/src/main/java/org/rcsb/project1/protein_chains_All_20150629_002251.seq";
//		 String outputfilename = "output.txt";
		 
		// Spark setup
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(Demo.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
		
		// start time
//		long start = System.nanoTime();
		
		// map sequence file to pairs (id, points)
		List<Tuple2<String, Point3d[]>> chains = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.0003, 123)
				.mapToPair(new HadoopToSimpleChainMapper()) // convert input to <pdbId.chainId, SimplePolymerChain> pairs
				.filter(t -> t._2.isProtein())
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(new LengthFilter(50,500)) // keep protein chains with 50 - 500 residues
				.collect();		
		
		// instantiations
		ArrayList<Tuple3<String, String, Point3d[]>> lib = new ArrayList<>();
		ArrayList<ArrayList<Tuple2<Double, Integer>>> comparisons = new ArrayList<>();
		ArrayList<Tuple2<Double, Integer>> templist = new ArrayList<>();
		SuperPositionQCP qcp = new SuperPositionQCP(true);
		ArrayList<Integer> skiplist = new ArrayList<>();
		ArrayList<Integer> numnulls = new ArrayList<Integer>();
		// int numnulls = 0;
		// int length = 8;
		int threshold = 1;
		
		// integer for comparison's index
		int ind = 0;
		
		// boolean for when to start adding to skiplist
		boolean go = false;
		
		// boolean for whether or not to add a fragment to the library
		boolean bool = true;
		





		for (Tuple2<String, Point3d[]> t: chains){
			for (int star = 0; star < t._2.length-length; star++){								// for each fragment
				Point3d[] fragment = new Point3d[length];
				for (int i = 0; i < length; i++) {
					fragment[i] = new Point3d(t._2[star+i]);
				}
				SuperPositionQCP.center(fragment);			
				
				Tuple3<String, String, Point3d[]> tup = 
						new Tuple3<String, String, Point3d[]>
							(t._1 + "." + star,
							lengthy(fragment),
							fragment);
				
				int u = lib.size();
				
				ind = 0;
				
				for (int i = 0; i < lib.size(); i++) {											// for each new fragment, loop once through the library
					if (!lib.isEmpty()) {
						if (lib.get(i)._2() == tup._2() && !skiplist.contains(i)) {
							// get (c)RMSD
							qcp.set(lib.get(i)._3(), tup._3());
							double q = qcp.getRmsd();
							
							Tuple2<Double, Integer> temptup = new Tuple2<Double, Integer>(q, ind);
							ind++;
							
							if (q < threshold) {
								bool = false;
								templist.clear();
								break;
							}
							templist.add(temptup);
							
							u -= numnulls.get(i);



							for (int a = 0; a < u; a++) {
								if (go) {
									if (comparisons.get(a).get(i)._1 - q > threshold) {
										go = false;
										break;
									}
									else if (!skiplist.contains(a)) {
										skiplist.add(a);
									}
								}
								else {
									if (q - comparisons.get(a).get(i)._1 <= threshold) {
										a -= 1;
										while (q - comparisons.get(a).get(i)._1 <= threshold) {
											a -= 1;
										}
										go = true;
									}
									if (q - comparisons.get((int) ((u-a)/2)).get(i)._1 > threshold) {
										a = (int) ((u-i)/2); // + 1 should be added automatically by for loop
									}
									else if (comparisons.get((int) ((u-a)/2)).get(i)._1 - q < threshold) {
										u = (int) ((u-a)/2);   															// CHECKCHECKCHECKCHECKCHECK
									}
								}
							}
						}
						else {
							templist.add(null);
						}
					}
				}
				skiplist.clear();




				if (bool == true) {
					starter: for (int vert = 0; vert < lib.size(); vert++) {
						for (int hor = vert; hor < lib.size(); hor++) {
							if (comparisons != null) {
								if (templist.get(vert) != null) {
									if (templist.get(vert)._1 <= comparisons.get(hor).get(vert)._1) {
										comparisons.get(hor).add(templist.get(vert));
									}
									else if (hor == lib.size()-1) {
										comparisons.get(hor).add(templist.get(vert));
									}
								}
								else {
									// templist(vert) is null
									comparisons.get(hor).add(templist.get(vert));
									numnulls.set(vert, numnulls.get(vert) + 1);					// CHECKCHECKCHECKCHECK
								}
							}
							else {
								//comparisons is completely null
								comparisons.add(templist);
								break starter;
							}
						}
					}
					templist.clear();
					lib.add(tup);
					//System.out.println("[" + tup._2() + "] - " + (lib.size()-1) + ": " + Arrays.toString(tup._3()));
				}
				bool = true;
			}
		}
		sc.close();




		PrintWriter writer = new PrintWriter("library.tsv", "UTF-8");
		for (Tuple3<String, String, Point3d[]> l: lib) {
			writer.print("PDBID.ChainID\tlength");
			for (int i=0; i<length; i++) {
				writer.print("\tPoint " + i);
			}
			writer.println();
			writer.print(l._1() + "\t [" + l._2() + "]");
			for (Point3d p: l._3()) {
				writer.print("\t" + p);
			}
			writer.println();
		}
		writer.close();
		System.out.println(lib.size());
		//System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}




	// FIND OPTIMAL LENGTH
	public static String lengthy(Point3d[] p){
		int round = 2;
		double i = Math.abs(p[p.length-1].distance(p[0]));
		i /= round;
		int base = (int) i * round;
		int top = ((int) i + 1) * round;
		return Integer.toString(base) + " - " + Integer.toString(top);
	}
}