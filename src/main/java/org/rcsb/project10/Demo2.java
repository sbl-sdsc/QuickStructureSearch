package org.rcsb.project10;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;

import javax.vecmath.Point3d;

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
public class Demo2 {

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
	    
	    JavaPairRDD<String, Point3d[]> mapToPair = sc
	    		.sequenceFile(chainFile, Text.class, WritableSegment.class) // read file with chains
	    		.mapToPair(t -> new Tuple2<String, WritableSegment> (new String(t._1.toString()), new WritableSegment(t._2)) ) // make a copy of the data
	    		.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()));


	  
	    sc.stop();
	    sc.close();
	    
		return 0;
	}
}
