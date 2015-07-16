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
	
	// make the data structure to store unique fragments
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException
	{
		 // arguments
		 String path = "/Users/grantsummers/Desktop/School/Internship/main/QuickStructureSearch/src/main/java/org/rcsb/project1/protein_chains_All_20150629_002251.seq";
//		 String outputfilename = "output.txt";
		 int length = 8;
		 
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
		
		
		
		
		
		
		
		ArrayList<Tuple3<String, String, Point3d[]>> lib = new ArrayList<>();
		
		// boolean for whether or not to add a fragment to the library
		boolean bool = true;
		
		// instantiate list of lists of RMSDs
		ArrayList<ArrayList<Tuple2<Double, Integer>>> comparisons = new ArrayList<>();
		
		// instantiate temp list
		ArrayList<Tuple2<Double, Integer>> templist = new ArrayList<>();
		
		// instantiate the SPQCP object
		SuperPositionQCP qcp = new SuperPositionQCP(true);
		
		// instantiate a skiplist
		ArrayList<Integer> skiplist = new ArrayList<>();
		
		// integer for comparison's index
		int ind = 0;
		
		bool go = false;
		
		// tuple2 to add to templis
		Tuple2<Double, Integer> temptup = new Tuple2<Double, Integer>(null, null); // if nulls don't work, use 0's
				
				
				
				
				
		
		
		
		
		for (Tuple2<String, Point3d[]> t: chains){
			for (int star = 0; star < t._2.length-length; star++){		
				
				// center fragment
				Point3d[] fragment = new Point3d[length];
				for (int i = 0; i < length; i++) {
					fragment[i] = new Point3d(t._2[star+i]);
				}
				SuperPositionQCP.center(fragment);			
				
				
				// create tuple3 fragment
				Tuple3<String, String, Point3d[]> tup = 
						new Tuple3<String, String, Point3d[]>
							(t._1 + "." + star,
							lengthy(fragment),
							fragment);
				int u = lib.size();
				if (!lib.isEmpty()) {
						if (lib.get(i)._2() == tup._2() && !skiplist.contains(i)) {
							
							// get (c)RMSD
							qcp.set(lib.get(i)._3(), tup._3());
							double q = qcp.getRmsd();
							
							temptup._1 = q;
							temptup._2 = ind;
							
/* CHECK THIS???*/							ind++;   //CHECK THIS?????
							
							// check if (c)RMSD is less than the threshold (1)
							if (q < 1) {
								bool = false;
								// if this tup isn't in lib, then we don't want to put useless info in comparisons
								templist.clear();
								break;
							}
							
							// Add RMSD to templist
							templist.add(temptup);
							
							System.out.println("before loop");
							
							// for (int k = i; k < lib.size(); k++) {
							
							
							// check which fragments in lib to skip w/binary search
							for (int a = int i; a < u; a++) {
								if (go) {
									if (comparisons.get(a).get(i) - q > 1) {
										//stop
										go = false;
										break;
									}
									else {
										// add index to skiplist if it's not already there
										if (!skiplist.contains(a)) {
											skiplist.add(a);
										}
									}
								}
								else {
									if (q - comparisons.get(a).get(i) <= 1) {
										//check behind (recursively), then start
										a -= 1;
										while (q - comparisons.get(a).get(i) <=1 ) {
											a -= 1;
										}
										go = true;
									}
									if (q - comparisons.get((int) ((u-a)/2)).get(i) > 1) {
										i = (int) ((u-i)/2); // + 1 should be added automatically by for loop
									}
									else if(comparisons.get((int) ((u-a)/2)).get(i) - q < 1){
										u = comparisons.get((int) ((u-a)/2)).get(i);
									}
								}
							}
								
								// if (Math.abs(q-comparisons.get(k).get(i)._1()) > 1) {  // check syntax w/ the ()
								// 	System.out.println("!!!added something to skiplist!!!");
								// 	// adds indices of fragments to skip to skiplist
								// 	skiplist.add(comparisons.get(k).get(i)._2);
								// }
								// else {
								// 	System.out.println"not dissimilar enough to add to skiplist");
								// }
							// }
						}
						else{templist.add(null);}
					}
					skiplist.clear();
				








					if (bool == true) {
						// inserts each item of templist into compariosons at the proper location (Increasing order in each row) w/binary search
						// note: null values are at the beginning
						starter: for (int vert = 0; vert < lib.size(); vert++) {
							for (int hor = vert; hor < lib.size(); hor++) {
								if (comparisons != null) {
									
									//compares elements of comparisons against elements of templist (horizontally)
									if (templist.get(vert) != null) {
										if (templist.get(vert) <= comparisons.get(hor).get(vert)) {
											// set comparisons(hor, vert) to templist(vert)
											comparisons.get(hor).add(templist(vert));
										}
										else if (hor == lib.size()-1) {
											comparisons.add(templist);
										}
									}
									else {
										// stuff goes here
									}
								}
								else {
									//comparisons is null
									comparisons.add(templist);
									break starter;
								}
							}
						}
						templist.clear();
						lib.add(tup);
	//					System.out.println("[" + tup._2() + "] - " + (lib.size()-1) + ": " + Arrays.toString(tup._3()));
					}
					
					
					
					
					
					
					
					
					
					
					bool = true;
				}
			}
		}
		sc.close();
		
		
		
		
		
		
		
		
		
		
		// Write the lib list to a text or csv file
//		PrintWriter writer = new PrintWriter("library.tsv", "UTF-8");
//		for(Tuple3<String, String, Point3d[]> l: lib){
//			writer.println(l._1() + "\t [" + l._2() + "]\t " + Arrays.toString(l._3()));
//		}
		
		// alternate write-out
		PrintWriter writer = new PrintWriter("library.tsv", "UTF-8");
		for(Tuple3<String, String, Point3d[]> l: lib){
			writer.print("PDBID.ChainID\tlength");
			for(int i=0; i<length; i++){
				writer.print("\tPoint " + i);
			}
			writer.println();
			writer.print(l._1() + "\t [" + l._2() + "]");
			for(Point3d p: l._3()){
				writer.print("\t" + p);
			}
			writer.println();
		}
		
		writer.close();
		System.out.println(lib.size());
		
		// prints time
//		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
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
