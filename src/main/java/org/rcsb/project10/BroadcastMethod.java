package org.rcsb.project10;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;

import org.apache.hadoop.io.Text;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.mmtf.spark.utils.SparkUtils;

import scala.Tuple2;
/**
 * This class performs pairwise structural alignments using broadcasting to eliminate
 * shuffles between servers.
 * 
 * @author Peter Rose
 *
 */
public class BroadcastMethod {

	public static void main(String[] args) {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		
		System.out.println("Broadcast Method");
		
		String chainFile = args[0];
		System.out.println("Chain file    : " + chainFile);
		
		String benchmarkFile = args[1] + File.separatorChar + "benchmark_" + timeStamp + ".csv";
	    System.out.println("Benchmark file: " + benchmarkFile);
		
	    long start = System.nanoTime();

		double fraction = 0.0002;
		int randomSeed = 1;
	    int combinations = run(chainFile, benchmarkFile, fraction, randomSeed);
	    
	    long end = System.nanoTime();
         
	    double time = (end-start)/1E9;
	    System.out.println("Time          : " + (end-start)/1E9 + " seconds"); 
	    System.out.println("Chain pairs   : " + combinations);
	    System.out.println("Throughput    : " + combinations/time + " pairs/seconds");
	}

	private static int run(String chainFile, String benchmarkFile,
			double fraction, int randomSeed) {
		JavaSparkContext sc = SparkUtils.getSparkContext();
	    
	    List<Tuple2<String, WritableSegment>> segments = sc
	    		.sequenceFile(chainFile, Text.class, WritableSegment.class) // read file with chains
	    		.mapToPair(t -> new Tuple2<String, WritableSegment> (new String(t._1.toString()), new WritableSegment(t._2)) ) // make a copy of the data
	    		.sample(false, fraction, randomSeed) // take a sample
	    		.collect();
	    
	    int combinations = (segments.size()*segments.size()-1)/2;
	    
	    // broadcast segment data to all nodes
	    final Broadcast<List<Tuple2<String,WritableSegment>>> data = sc.broadcast(segments); 

	    // get indices for all unique pairs given the number of segments
	    JavaPairRDD<Integer, Integer> pairs = SparkUtils.getComparisonMatrix(segments.size()); 
	    
	    // run structural alignment and save results in .csv file

	    pairs
	    	.repartition(8)
	    	.map(new StructuralAlignmentMapper(data))
	    	.saveAsTextFile(benchmarkFile);
	 
	    sc.stop();
	    sc.close();
	    
		return combinations;
	}
}
