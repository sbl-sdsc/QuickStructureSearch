package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class Demo1 {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]").setAppName(Demo1.class.getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);

		long start = System.nanoTime();
		
		// read sequence file and map sequence length to an RDD
		List<Tuple2<String, Point3d[]>> list = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.002)
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate array> pairs
				.filter(new GapFilter(0, 0)) // filter chains with zero gap length and zero gaps
				.collect();
		
//		for (Tuple2<String, Point3d[]> t: list) {
//			System.out.println(t._1 + ": " + Arrays.toString(t._2));
//		}

		for (int i = 0; i < list.size()-1; i++) {
			Point3d[] chain1 = list.get(i)._2;
			Point3d[] chain2 = list.get(i+1)._2;

			if (chain1.length == chain2.length) {
				String chainId1 = list.get(i)._1;
				String chainId2 = list.get(i+1)._1;
	
				
				SuperPositionQCP qcp = new SuperPositionQCP();
				qcp.set(chain1, chain2);
				double rmsd = qcp.getRmsd();
				
				System.out.println("Rmsd of :" + chainId1 + " - " + chainId2 + " = " + rmsd) ;
				
			}
		}

		

		
		
		
		
		
		sc.close();
		
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
}
