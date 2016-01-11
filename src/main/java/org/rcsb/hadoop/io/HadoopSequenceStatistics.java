package org.rcsb.hadoop.io;

import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.Comparator;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file with
 * PDBId.ChainId/SimplePolymerChain pairs. This demo program calculates the size of each
 * SimplePolymerChain data structure and calculates basic statistics.
 * 
 * @author  Peter Rose
 */
public class HadoopSequenceStatistics implements Serializable {    
	private static final long serialVersionUID = 1L;
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static final Logger logger = LoggerFactory.getLogger(HadoopSequenceStatistics.class);
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		// Hadoop sequence input file name
		String path = args[0]; 

		// Setup Spark
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(HadoopSequenceStatistics.class.getSimpleName())
				.set("spark.ui.showConsoleProgress", "false")
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		
		// Setup custom serialization for classes used in this application
		Class<?>[] classes = {SimplePolymerChain.class};
		conf.registerKryoClasses(classes);		
		
		JavaSparkContext sc = new JavaSparkContext(conf);

		long start = System.nanoTime();
		
		// read sequence file and size (in bytes) of the underlying SimplePolymerChain data structure
		JavaRDD<Integer> len = sc
				.sequenceFile(path, Text.class, SimplePolymerChain.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
			    .map(t -> t._2.size())
				.cache(); // cache since we are using the JavaRDD multiple times below
		
		// Calculate some basic statistics on the size of the data structure
		long count = len.count();
		
		logger.info(count + " records processed");
		
	    long min = len.min(Comparator.naturalOrder());
		long max= len.max(Comparator.naturalOrder());
		long total = len.reduce((a,b) -> a + b);
		long mean = total/count;
		
		long end = System.nanoTime();
		
		// shutdown Spark
		sc.stop();
		sc.close();
				
        System.out.println("Total chains         : " + count);
        System.out.println("Size in bytes ");
		System.out.println("Total  : " + total);	
		System.out.println("Mean   : " + mean);
		System.out.println("Minimum: " + min);
		System.out.println("Maximum: " + max);
		
		System.out.println("Time: " + (end - start)/1E9 + " sec.");
	}
}
