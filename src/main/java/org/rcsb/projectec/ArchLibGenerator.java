package org.rcsb.projectec;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.Accumulable;
import org.apache.spark.SparkConf;
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
		
		System.out.println("ArchLibGen");
		
		String chainFile = args[0];
		System.out.println("Chain file    : " + chainFile);
		
	    long start = System.nanoTime();

	//	double fraction = 0.0002;
	//	int randomSeed = 1;
	 //   int combinations = run(chainFile, benchmarkFile, fraction, randomSeed);
	    getAllFragments(chainFile,.01,1);
	    long end = System.nanoTime();
         
	    double time = (end-start)/1E9;
	    System.out.println("Time          : " + (end-start)/1E9 + " seconds"); 
	//    System.out.println("Chain pairs   : " + combinations);
	  //  System.out.println("Throughput    : " + combinations/time + " pairs/seconds");
	}

	private static JavaRDD<Point3d[]> getAllFragments(String chainFile,
			double fraction, int randomSeed){
		// setup spark

		SparkConf conf = new SparkConf()

		.setMaster("local[*]")

		.setAppName("ArchLibGenerator")

		.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")

		.set("spark.driver.maxResultSize", "3500m");


		JavaSparkContext sc = new JavaSparkContext(conf);
	    
	    JavaPairRDD<String, Point3d[]> mapToPair = sc
	    		.sequenceFile(chainFile, Text.class, WritableSegment.class) // read file with chains
	    		.mapToPair(t -> new Tuple2<String, WritableSegment> (new String(t._1.toString()), new WritableSegment(t._2)) ) // make a copy of the data
	    		.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
	    		.sample(false, fraction, randomSeed);
	
	    JavaRDD<Point3d[]>fragments = mapToPair.map(t->t._2()).flatMap(new FlatMapToFragments(6));
//	    fragments.foreach(t-> System.out.println(Arrays.toString(t)));

		//issue with adding, maybe just use list... make for loop and filter out stuff
		List<Point3d[]> prelib = new ArrayList<>();
		List<Point3d[]> lib = Collections.synchronizedList(prelib);
		Accumulable<List<Point3d[]>,Point3d[]> accLib = new Accumulable<>(lib,new AccumuableList());
		
		fragments.foreach(t->accLib.add(t));
		List<Point3d[]> pointList = accLib.value();
		
		System.out.println(pointList.size());
	
		


	    sc.stop();
	    sc.close();
//	    
		return fragments;
	}
	
//	private static List<Point3d[]> createLibrary(JavaRDD<Point3d[]> allFragments){
//		
//		List<Point3d[]> lib = new ArrayList<>();
////		for (Point3d[] fragment : allFragments) {
////			
////		}
//		
//		Accumulable<List<Point3d[]>,Point3d[]> accLib = new Accumulable<>(lib,new AccumuableList());
//		
//		
//		return null;
//		
//	}
	
	
}
