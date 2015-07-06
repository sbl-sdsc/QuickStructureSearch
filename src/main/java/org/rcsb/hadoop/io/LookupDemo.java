package org.rcsb.hadoop.io;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;

import scala.Tuple2;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class LookupDemo {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(LookupDemo.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		
		JavaSparkContext sc = new JavaSparkContext(conf);
		

		long start = System.nanoTime();
		
		// if you need both the coordinates and the sequences, use this section of code
		// read sequence file and map to PdbId.chainId, SimplePolymerChain pairs
        JavaPairRDD<String, SimplePolymerChain> chains = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
//				.sample(false, 0.01, 123)
				.mapToPair(new HadoopToSimpleChainMapper()) // convert input to <pdbId.chainId, SimplePolymerChain> pairs
				.filter(t -> t._2.isProtein())
				.cache();
        long start1 = System.nanoTime();
        int n = 0;
        
        List<String> sample = chains.sample(false,  0.01, 123).map(t -> t._1).collect();
        for (String s: sample) {
        	List<SimplePolymerChain> lookup = chains.lookup(s);
        	for (SimplePolymerChain c: lookup) {
        		n++;
        	}
        }
        long end = System.nanoTime();
        System.out.println("Looked up: " + n);
        System.out.println("Time: " + (end-start1)/1E6 + " ms");
        System.out.println("Time per entry: " + (end-start1)/1E6/n + " ms");
		
		// if you need just the coordinates, use this section of code
		// read sequence file and map to PdbId.chainId, C-alpha coordinate pairs
//		List<Tuple2<String, Point3d[]>> coordinates = sc
//				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
//				.sample(false, 0.01, 123)
//				.mapToPair(new HadoopToSimpleChainMapper()) // convert input to <pdbId.chainId, protein sequence> pairs
//				.filter(t -> t._2.isProtein())
//				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
//				.collect();
//
//		for (Tuple2<String, Point3d[]> t: coordinates) {
//			System.out.println(t._1 + ": " + Arrays.toString(t._2));
//		}

		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
}
