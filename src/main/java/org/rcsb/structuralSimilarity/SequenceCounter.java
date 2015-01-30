package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.util.Comparator;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class SequenceCounter {    
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 2; // Spark recommends 2-3 tasks per thread
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]").setAppName(SequenceCounter.class.getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);

//		JavaPairRDD<Text, ArrayWritable> seq = sc.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD);      
//		JavaRDD<Long> len = seq.map(s -> new Long(((IntWritable)s._2.get()[0]).get()));
//		len.cache();

		long start = System.nanoTime();
		
		// compact version
		JavaRDD<Long> len = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.map(s -> new Long(((IntWritable)s._2.get()[0]).get()))
				.cache(); // cache since we are using the JavaRDD multiple times below
		
		long chainCount = len.count();
		long minChainLength = len.min(Comparator.naturalOrder());
		long maxChainLength = len.max(Comparator.naturalOrder());
		
		long residueCount = len.reduce((a,b) -> a + b);
		
		sc.stop();
		sc.close();
				
        System.out.println("Total chains         : " + chainCount);
		System.out.println("Total residues       : " + residueCount);	
		System.out.println("Average chain length : " + residueCount/chainCount);
		System.out.println("Minimum chain length : " + minChainLength);
		System.out.println("Maximum chain length : " + maxChainLength);
		
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
}
