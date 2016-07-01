package org.rcsb.projectec;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.mmtf.spark.utils.SparkUtils;
import org.rcsb.project10.StructuralAlignmentMapper;

import org.rcsb.project10.WritableSegment;
import org.rcsb.structuralAlignment.SuperPositionQCP; // yay

import scala.Tuple2;
/**
 * This class performs pairwise structural alignments using broadcasting to eliminate
 * shuffles between servers.
 * 
 * @author Peter Rose
 *
 */
public class ArchLibGenerator {
	public static void main(String[] args) {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		
		System.out.println("Broadcast Method");
		
		String chainFile = args[0];
		System.out.println("Chain file    : " + chainFile);
		
	    long start = System.nanoTime();

	//	double fraction = 0.0002;
	//	int randomSeed = 1;
	 //   int combinations = run(chainFile, benchmarkFile, fraction, randomSeed);
	    getAllFragments(chainFile);
	    long end = System.nanoTime();
         
	    double time = (end-start)/1E9;
	    System.out.println("Time          : " + (end-start)/1E9 + " seconds"); 
	//    System.out.println("Chain pairs   : " + combinations);
	  //  System.out.println("Throughput    : " + combinations/time + " pairs/seconds");
	}

	private static JavaRDD<Point3d[]> getAllFragments(String chainFile){
		JavaSparkContext sc = SparkUtils.getSparkContext();
	    
	    JavaPairRDD<String, Point3d[]> mapToPair = sc
	    		.sequenceFile(chainFile, Text.class, WritableSegment.class) // read file with chains
	    		.mapToPair(t -> new Tuple2<String, WritableSegment> (new String(t._1.toString()), new WritableSegment(t._2)) ) // make a copy of the data
	    		.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()));
	
	    JavaRDD<Point3d[]>fragments = mapToPair.map(t->t._2()).flatMap(new FlatMapToFragments(6));
	    fragments.foreach(t-> System.out.println(Arrays.toString(t)));

	    
	  
	    sc.stop();
	    sc.close();
	    
		return fragments;
	}
	
	
}
