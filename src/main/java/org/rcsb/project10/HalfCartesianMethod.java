package org.rcsb.project10;

import java.util.Arrays;

import org.apache.hadoop.io.Text;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.mmtf.spark.utils.SparkUtils;

import scala.Tuple2;

public class HalfCartesianMethod {

	public static void main(String[] args) {
	    long start = System.nanoTime();

	    JavaSparkContext sc = SparkUtils.getSparkContext();
	    
	    JavaPairRDD<String, WritableSegment> segments = sc
	    		.sequenceFile("/Users/peter/Data/MMTF/x-rayChains.seq", Text.class, WritableSegment.class)
	    		.mapToPair(t -> new Tuple2<String, WritableSegment> (new String(t._1.toString()), new WritableSegment(t._2)) )
	            .sample(false, 0.0001, 1);
//	    System.out.println("segments:  " + segments.count());
	 
	    JavaPairRDD<Tuple2<String, WritableSegment>, Tuple2<String, WritableSegment>> halfCartesian = SparkUtils.getHalfCartesian(segments);
	    JavaPairRDD<String, Float[]> scores = halfCartesian.mapToPair(new SegmentPairMapperToTmMapper()).cache();
	    scores.foreach(t -> System.out.println(t._1 + " " + Arrays.toString(t._2)));
	    System.out.println("pairs:  " + scores.count());
	    
//	    JavaRDD<Integer> scores = halfCartesian.map(t -> Math.abs(t._1._1.length() - t._2._1.length()));
//	    scores.foreach(t -> System.out.println(t + " " + Arrays.toString));
	    
	    long end = System.nanoTime();        
	    System.out.println("Time: " + (end-start)/1E9);
	    
	    sc.close();
	}
}
