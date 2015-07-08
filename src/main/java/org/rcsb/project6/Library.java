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
		
		List<Tuple3<String, String, Point3d[]>> lib = new ArrayList<>();
		
		// boolean for whether or not to add a fragment to the library
		boolean bool = true;
		
		// instantiate the SPQCP object
		SuperPositionQCP qcp = new SuperPositionQCP(true);
		
		for (Tuple2<String, Point3d[]> t: chains) {
			for(int star=0; star<t._2.length-length; star++)
			{				
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
					check: for(Tuple3<String, String, Point3d[]> l: lib){
						if(l._2() == tup._2()){
							qcp.set(l._3(), tup._3());
							double q = qcp.getRmsd();
							if(q<1){
								bool = false;
								break check;
							}
						}
					}
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
		PrintWriter writer = new PrintWriter("library.csv", "UTF-8");
		for(Tuple3<String, String, Point3d[]> l: lib){
			writer.print(l._1() + ", [" + l._2() + "]");
			for(Point3d p: l._3()){
				writer.print(", \"" + p + "\"");
			}
			writer.println();
		}
		
		writer.close();
		
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
