package org.rcsb.project10;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.mmtf.spark.utils.SparkUtils;

import scala.Tuple2;

public class SegmentReader {

	public static void main(String[] args) {
	    long start = System.nanoTime();

	    SparkConf conf = new SparkConf().setMaster("local[*]").setAppName(SparkUtils.class.getSimpleName()); 
	    JavaSparkContext sc = new JavaSparkContext(conf);
	    
	    JavaPairRDD<String, WritableSegment> segments = sc
	    		.sequenceFile("/Users/peter/Data/MMTF/x-rayChains.seq", Text.class, WritableSegment.class)
	    		.mapToPair(t -> new Tuple2<String, WritableSegment> (new String(t._1.toString()), t._2) );
	    
	    segments = segments.filter(t -> t._2.getSequence().contains("HHHHH"));	    
	    System.out.println("Count: " + segments.count());
	    
//	    structure.foreach(t -> System.out.println(t._1));
//	    structure.foreach(t -> System.out.println(Arrays.toString(t._2.getCoordinates())));
	 
	    long end = System.nanoTime();
         
	    System.out.println("Time: " + (end-start)/1E9);
	    sc.close();
	}
}
