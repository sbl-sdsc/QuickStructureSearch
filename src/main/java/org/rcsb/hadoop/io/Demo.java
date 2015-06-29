package org.rcsb.hadoop.io;

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
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.SeqToChainMapper;

import scala.Tuple2;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class Demo {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]").setAppName(Demo.class.getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);

		long start = System.nanoTime();

		// read sequence file and map to PdbId.chainId, C-alpha coordinate pairs
		List<Tuple2<String, Point3d[]>> coordinates = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.001, 123)
				.mapToPair(new HadoopToCoordinateMapper()) // convert input to <pdbId.chainId, CA coordinate array> pairs
				.collect();

		for (Tuple2<String, Point3d[]> t: coordinates) {
			System.out.println(t._1 + ": " + Arrays.toString(t._2));
		}

		// read sequence file and map to PdbId.chainId, protein sequence pairs
		List<Tuple2<String, String>> sequences = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.001, 123)
				.mapToPair(new HadoopToSequenceMapper()) // convert input to <pdbId.chainId, protein sequence> pairs
				.collect();

		for (Tuple2<String, String> t: sequences) {
			System.out.println(t._1 + ": " + t._2);
		}

		// read sequence file and map to PdbId.chainId, SimplePolymerChain pairs
		List<Tuple2<String, SimplePolymerChain>> chains = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.1, 123)
				.mapToPair(new HadoopToSimpleChainMapper()) // convert input to <pdbId.chainId, protein sequence> pairs
				.filter(t -> t._2.isProtein())
				.collect();

		for (Tuple2<String, SimplePolymerChain> t: chains) {
			System.out.println(t._1 + ": " + t._2);
		}
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
}
