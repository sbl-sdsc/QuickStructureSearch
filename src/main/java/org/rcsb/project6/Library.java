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
		ArrayList<ArrayList<Double>> comparisons = new ArrayList<ArrayList<Double>>();
		
		// instantiate temp list
		ArrayList<Double> templist = new ArrayList<Double>();
		
		// instantiate the SPQCP object
		SuperPositionQCP qcp = new SuperPositionQCP(true);
		
		// instantiate a skiplist
		ArrayList<Integer> skiplist = new ArrayList<Integer>();
		
		// instantiate a sort list
		ArrayList<Double> sortlist = new ArrayList<Double>();
		
		// pseudocode
		/*
		 * instantiations
		 * 
		 * for each chain (label each one t)
		 *     for each fragment in t
		 *         center all points in t about centroid (set centroid = (0, 0, 0))
		 *         create a 3-tuple of the centered fragment
		 *         if library is not empty:
		 *             (libcheck) for each fragment in the library (each has index i):
		 *                 if length of lib[i] is approximately the length of the 3-tuple /////////AND i is NOT in skiplist:
		 *                     let q = the rmsd of the 3-tuple and lib[i]
		 *                     for all rmsds in the row of lib[i], from i on, in comparisons:
		 *                         if absolute value(q - rmsd(in compositions)>1):
		 *                             add index of other fragment in compositions to skiplist
		 *                     add q to templist
		 *                     if q was less than 1:
		 *                         bool = false
		 *                         clear templist (b/c we don't want to add lists of values that correspond to a point that was discarded to comparisons
		 *                         end the libcheck for loop
		 *                 else (lengths not the same):
		 *                     add a null value to templist
		 *             add templist (now populated with a bunch of rmsd values) to comparisons
		 *             clear templist to prep for next loop
		 *             clear the skiplist
		 *         if bool = true:
		 *             add the 3-tuple to lib
		 *         set bool = true for the next fragment check (run another 3-tuple against the library)
		*/
		
		for (Tuple2<String, Point3d[]> t: chains){
			for(int star=0; star<t._2.length-length; star++){				
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
				
				if(!lib.isEmpty()){
					check: for(int i = 0; i<lib.size(); i++){
						if(lib.get(i)._2() == tup._2() && !skiplist.contains(i)){
							
							// get (c)RMSD
							qcp.set(lib.get(i)._3(), tup._3());
							double q = qcp.getRmsd();
							
							// Add RMSD to templist
							templist.add(q);
							
							// check if (c)RMSD is less than the threshold (1)
							if(q<1){
								bool = false;
								// if this tup isn't in lib, then we don't want to put useless info in comparisons
								templist.clear();
								break check;
							}
							
							System.out.println("before loop");
							
							// check which fragments in lib to skip
							for(int k = i; k<lib.size();k++){
								System.out.println("k loop");
								if(comparisons.get(k).get(i)!=null){
									System.out.println("comparisons is not null");
									if(Math.abs(q-comparisons.get(k).get(i))>1){
										System.out.println("made it!");
										
										// adds indices of fragments to skip to skiplist
										skiplist.add(k);
									}
								}
							}
						}
						else{templist.add(null);}
					}
					comparisons.add(templist);
					templist.clear();
					skiplist.clear();
				}
				
				if(bool == true){
					lib.add(tup);
//					System.out.println("[" + tup._2() + "] - " + (lib.size()-1) + ": " + Arrays.toString(tup._3()));
				}
				bool = true;
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
	
	public static String lengthy(Point3d[] p){
		int round = 2;
		double i = Math.abs(p[p.length-1].distance(p[0]));
		i /= round;
		int base = (int) i * round;
		int top = ((int) i + 1) * round;
		return Integer.toString(base) + " - " + Integer.toString(top);
	}
}
