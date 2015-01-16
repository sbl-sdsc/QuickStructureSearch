package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.util.Comparator;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;

/**
 * 
 * @author  Peter Rose
 */
public class SequenceCounter 
{    
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = "/Users/peter/Data/testdata/test.seq";
		
		long start = System.nanoTime();

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		int cores = 2;
		SparkConf conf = new SparkConf().setMaster("local[" + cores + "]").setAppName(SequenceCounter.class.getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);

		JavaPairRDD<Text, ArrayWritable> seq = sc.sequenceFile(path, Text.class, ArrayWritable.class,cores);      
	
		JavaRDD<Long> len = seq.map(s -> new Long(s._2.get().length/3));
		len.cache();
		
		long residueCount = len.reduce((a,b) -> a + b);
		
		long chainCount = len.count();
		long minChainLength = len.min(Comparator.naturalOrder());
		long maxChainLength = len.max(Comparator.naturalOrder());
		sc.stop();
		sc.close();
				
        System.out.println("Total chains: " + chainCount);
		System.out.println("Total residues count : " + residueCount);	
		System.out.println("Average chain length : " + residueCount/chainCount);
		System.out.println("Minimum chain length : " + minChainLength);
		System.out.println("Maximum chain length : " + maxChainLength);
		
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
}
