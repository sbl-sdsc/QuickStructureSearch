package org.rcsb.project10;

import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.SequenceFileOutputFormat;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.mmtf.spark.utils.SparkUtils;

import scala.Tuple2;

public class Repartitioner {

	public static void main(String[] args) {
	    long start = System.nanoTime();

	    SparkConf conf = new SparkConf().setMaster("local[*]")
				.setAppName(SparkUtils.class.getSimpleName()); 
	    JavaSparkContext sc = new JavaSparkContext(conf);
	    
	    sc.sequenceFile("/Users/peter/Data/MMTF/reduced", Text.class, BytesWritable.class)
	    		.coalesce(12) 
                .mapToPair(t -> new Tuple2<Text, BytesWritable>(new Text(t._1.toString()), new BytesWritable(t._2.getBytes()) ))
                .saveAsHadoopFile("/Users/peter/Data/MMTF/reduced8", Text.class, BytesWritable.class,  SequenceFileOutputFormat.class);
	    
	    long end = System.nanoTime();
         
	    System.out.println("Time: " + (end-start)/1E9);
	    sc.close();
	}
}
