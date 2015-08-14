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

import org.rcsb.structuralAlignment.SuperPositionQCP;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;

import scala.Tuple2;
import scala.Tuple3;


/* 
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

	/**
	 * Main method
	 * @param  args                         no arguments necessary - all specified. Make sure to modify 
	 *                                      the path below to fit your computer!
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 */
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException
	{
		
		String path = "/Users/grantsummers/Desktop/School/Internship/main/QuickStructureSearch/src/main/java/org/rcsb/project1/protein_chains_All_20150629_002251.seq";

		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(Demo.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
		
		long start = System.nanoTime();

		List<Tuple2<String, Point3d[]>> chains = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.003, 123)
				.mapToPair(new HadoopToSimpleChainMapper()) 
				.filter(t -> t._2.isProtein())
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
				.filter(new GapFilter(0, 0)) 
				.filter(new LengthFilter(50,500)) 
				.collect();		
		
		ArrayList<Tuple3<String, String, Point3d[]>> lib = new ArrayList<>();
		SuperPositionQCP qcp = new SuperPositionQCP(true);
		int threshold = 1;
		boolean bool = true;
		int counter = 0;
		ArrayList<String> title = new ArrayList<>();
		ArrayList<Integer> freq = new ArrayList<>();
		int number;
		double q;
		ArrayList<ArrayList<Double>> comparisons = new ArrayList<>();
		ArrayList<Double> templist = new ArrayList<>();
		ArrayList<Integer> skiplist = new ArrayList<>();

		
		for (Tuple2<String, Point3d[]> t: chains){
			for (int star = 0; star < t._2.length-length; star++){
				counter ++;
				
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

				if (lib.size() > 0) {
					libloop: for (int i = 0; i < lib.size(); i++) {
						if (lib.get(i)._2().equals(tup._2()) && !skiplist.contains(i)) {
							qcp.set(lib.get(i)._3(), tup._3());
							q = qcp.getRmsd();
							for (int j = i; j < comparisons.size(); j++) {
								if (i == lib.size() - 1) {
									System.out.println("i and end");
									break;
								}
								else if (comparisons.get(j).size() > i && Math.abs(q - comparisons.get(j).get(i)) >= 1) {
									skiplist.add(j);
									System.out.println("something added to  skipplist!");
								}
							}
							if (q < threshold) {
								bool = false;
								templist.clear();
								int key = searcher(title, lib.get(i)._1());
								if (key != -1) {
									number = freq.get(key) + 1;
									freq.set(key, number);
								}
								else {
									title.add(tup._1());
									freq.add(1);
								}
								break libloop;
							}
							templist.add(q);
						}
						else {
							templist.add(null);
						}
					}
					skiplist.clear();


					if (bool == true) {
						title.add(tup._1());
						freq.add(1);
						lib.add(tup);
						comparisons.add(templist);
						templist.clear();
					}
					else {
						bool = true;
					}
				}
				else {
					lib.add(tup);
					freq.add(1);
					title.add(tup._1());
				}
			}
		}
		sc.close();

		PrintWriter writer = new PrintWriter("library.txt", "UTF-8");
		writer.print("PDBID.ChainID\tlength");
		for (int i=1; i<=length; i++) {
			writer.print("\tPoint " + i);
		}
		writer.print("\tFrequency");
		writer.println();
		for (int i=0; i < lib.size(); i++) {
			writer.print(lib.get(i)._1() + "\t [" + lib.get(i)._2() + "]");
			for (Point3d p: lib.get(i)._3()) {
				writer.print("\t" + p);
			}
			writer.print("\t" + freq.get(i));
			writer.println();
		}
		writer.close();
		System.out.println("library size: " + lib.size());
		System.out.println("initial size: " + counter);
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
	

	/**
	 * Finds the index of a string in an arraylist, to be used in conjunction with another arraylist.
	 * @param  a Arraylist to search through
	 * @param  b Search argument
	 * @return   index of b in a - if it exists, otherwise, returns -1
	 */
	public static Integer searcher (ArrayList<String> a, String b) {
		for (int s = 0; s < a.size(); s++) {
			if (a.get(s).equals(b)) {
				return s;
			}
		}
		return -1;
	}

	/**
	 * Returns a string in the form of "x - y", where x is the length of p rounded down
	 * to the nearest @round and y is the length of p rounded up to the nearest @round.
	 * @param p point 3d array that this method finds the length boundaries of.
	 * @return String with length of p rounded up and down to the nearest @round.
	 */
	public static String lengthy (Point3d[] p) {
		int round = 2;
		double i = Math.abs(p[p.length-1].distance(p[0]));
		i /= round;
		int base = (int) i * round;
		int top = ((int) i + 1) * round;
		return Integer.toString(base) + " - " + Integer.toString(top);
	}
}