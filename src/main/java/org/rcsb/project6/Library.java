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
//import org.rcsb.hadoop.io.SimplePolymerChain;
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
//		String path = "/Users/MellissaSummers/QuickStructureSearch/src/main/java/org/rcsb/project6/protein_chains_All_20150629_002251.seq";
//		String path = "~/QuickStructureSearch/src/main/java/org/rcsb/project6/protein_chains_All_20150629_00251.seq";
		
//		String outputfilename = "output.txt";
		 
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
		// lib is the library (a list) of unique fragments
		ArrayList<Tuple3<String, String, Point3d[]>> lib = new ArrayList<>();
		// 2d arraylist of tuple2's that contains RMSDs between fragments in lib
		// the second value in the tuple2 contains the horizontal index so that order can be
		// found even after they are horizontally sorted (each row individually sorted)
		ArrayList<ArrayList<Tuple2<Double, Integer>>> comparisons = new ArrayList<>();
		ArrayList<Tuple2<Double, Integer>> templist = new ArrayList<>();
		SuperPositionQCP qcp = new SuperPositionQCP(true);
		// skiplist is a list of fragments in library to skip comparisons of
		ArrayList<Integer> skiplist = new ArrayList<>();
		// numnulls is a list of integers, where each integer represents the number
		// of null values in comparisons in the corresponding column
		ArrayList<Integer> numnulls = new ArrayList<Integer>();
		// int length = 8;
		int threshold = 1;
		// integer for comparison's index
		int ind;
		// boolean for when to start adding to skiplist
		boolean go = false;
		// boolean for whether or not to add a fragment to the library
		boolean bool = true;
		int u;
		





		for (Tuple2<String, Point3d[]> t: chains){
			for (int star = 0; star < t._2.length-length; star++){								// for each fragment
				
				// code applied to (each) fragment with length 8 and start position star on chain t.
				
				Point3d[] fragment = new Point3d[length];
				for (int i = 0; i < length; i++) {
					fragment[i] = new Point3d(t._2[star+i]);
				}
				SuperPositionQCP.center(fragment);
				
				// fragment now has point3d[] = fragment, which is centered about centroid.
				
				Tuple3<String, String, Point3d[]> tup = 
						new Tuple3<String, String, Point3d[]>
							(t._1 + "." + star,
							lengthy(fragment),
							fragment);
				
				
				if (lib.size() > 0) {
					u = lib.size();
					ind = lib.size() - 1;
					
					// loops through lib to check if tup is unique relative to lib's fragments
					libloop: for (int i = 0; i < lib.size(); i++) {
						if (lib.get(i)._2().equals(tup._2()) && !skiplist.contains(i)) {
							qcp.set(lib.get(i)._3(), tup._3());
							double q = qcp.getRmsd();
							
							Tuple2<Double, Integer> temptup = new Tuple2<Double, Integer>(q, ind);
							ind--;
							System.out.println("temptup - [q, ind]: [" + q + ", " + ind + "]");


							
							if (q < threshold) {
								System.out.println("skipped due to similarity");
								bool = false;
								if(!templist.isEmpty()){
									templist.clear();
								}
								break libloop;
							}
							templist.add(temptup);
							
							if (numnulls.size() >= i && numnulls.get(i) != null) { 
								u -= numnulls.get(i);
							}
							
							
							
							
							 for (int a = 0; a < u; a++) {
//								System.out.println("comparisons: " + comparisons.size());
//								System.out.println("comparisons[a]: " + comparisons.get(a).size());
//								System.out.println("(int) ((u-1-a)/2): " + ((int) ((u-1-a)/2)));
//								System.out.println("comparisons[half]: " + comparisons.get((int) ((u-1-a)/2)).size());
//								System.out.println("u: " + u);
//								System.out.println("a: " + a);
//								System.out.println("i: " + i);
								if (comparisons.get(a).get(i) != null && comparisons.get((int) ((u-1-a)/2)).get(i) != null){
									System.out.println("difference (q - comparisons): " + (q-comparisons.get(a).get(i)._1));
									double difference = q - comparisons.get(a).get(i)._1;
									double halfdiff = q - comparisons.get((int) ((u-1-a)/2)).get(i)._1;
									if (go) {
										System.out.println("we are go");
										if (-difference > threshold) {
											System.out.println("end go");
											go = false;
											break;
										}
										else if (!skiplist.contains(a)) {
											skiplist.add(a);
											System.out.println("Skiplist: " + a);
										}
									}
									else {
										if (difference <= threshold) {
											System.out.println("start");
											a -= 1;
											while (difference <= threshold) {
												a -= 1;
											}
											go = true;
										}
										if (halfdiff > threshold) {
											System.out.println("1");
											a = (int) ((u-1-a)/2);
										}
										else if (-halfdiff > threshold) {
											System.out.println("2");
											u = (int) ((u-1-a)/2);
										}
										else if (halfdiff == threshold) {
											System.out.println("3 - start");
											a = (int) ((u-1-a)/2) - 1;
											while (difference <= threshold) {
												a -= 1;
											}
											go = true;
											
										}
									}
								}
								else {
									System.out.println("derp: failed triangle inequality (prolly a sort - insertion - problem)");
								}
							}
							System.out.println("I'm' outta the loop");
						}
						else {
							templist.add(null);
						}
					}
					skiplist.clear();

					System.out.println("templist: " + templist.size());
					if (bool == true) {
						vert: for (int vert = 0; vert < lib.size(); vert++) {
							hor: for (int hor = vert; hor < lib.size(); hor++) {
								if (comparisons != null && !comparisons.isEmpty()) {
									if (templist.get(vert) != null) {
										if (templist.get(vert)._1 <= comparisons.get(hor).get(vert)._1 || hor == lib.size() - 1) {
											comparisons.get(hor).add(templist.get(vert));
											break hor;
										}
									}
									else {
										//templist(vert) is null
										System.out.println("comparisons: " + comparisons.size());
										System.out.println("lib-1: " + (lib.size()-1));
										comparisons.get(lib.size()-1).add(vert, null);
										if (numnulls != null && numnulls.size() > vert) {
											numnulls.set(vert, numnulls.get(vert) + 1);
										}
										else if (numnulls.size() < vert - 1) {
//											while (numnulls.size() < vert - 2) {
//												numnulls.add(null);
//											}
											numnulls.add(1);
										}
									}
								}
								else {
									//comparisons is completely null
									comparisons = new ArrayList<ArrayList<Tuple2<Double, Integer>>>();
									comparisons.add(templist);
									break vert;
								}
							}
						}
						templist.clear();
						lib.add(tup);
					}
					else {
						bool = true;
					}
				}
				else {
					lib.add(tup);
				}
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
		System.out.println("library size: " + lib.size());
		
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