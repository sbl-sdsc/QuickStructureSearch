package org.rcsb.project10;

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
public class ArchLibGeneratorPR {
	public static void main(String[] args) {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		
		System.out.println("ArchLibGen");
		
		String chainFile = args[0];
		System.out.println("Chain file    : " + chainFile);
		
	    long start = System.nanoTime();

		int fragmentSize = 8;
		double rmsdThreshold = 2.0;

	    getAllFragments(chainFile, fragmentSize, rmsdThreshold);
	    
	    long end = System.nanoTime();
         
	    System.out.println("Time          : " + (end-start)/1E9 + " seconds"); 
	}

	private static void getAllFragments(String chainFile, int fragmentSize, double rmsdThreshold){
		
		// setup spark
		SparkConf conf = new SparkConf()
		.setMaster("local[*]")
		.setAppName("ArchLibGenerator")
		.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
	
		// read protein and cut into fragments (sliding window approach)
	    JavaRDD<Point3d[]> fragments = sc
	    		.sequenceFile(chainFile, Text.class, WritableSegment.class) // read file with chains
	    		.map(t -> t._2.getCoordinates()) // get the coordinates of the protein chains
	    		.repartition(1) // create a single partition to generate a single fragment library (this cannot be done in parallel!)
	    		.flatMap(new FlatMapToFragmentsEC(8)); // flatmap to fragments

//	    fragments.foreach(t-> System.out.println(Arrays.toString(t)));

	    // Create library of fragment archetypes
		List<Point3d[]> prelib = new ArrayList<>();
		List<Point3d[]> lib = Collections.synchronizedList(prelib);
		Accumulable<List<Point3d[]>,Point3d[]> accLib = new Accumulable<>(lib,new AccumuableListPR(rmsdThreshold));
		
		fragments.foreach(t->accLib.add(t));
		
		List<Point3d[]> archetypes = accLib.value();
		
		System.out.println(archetypes.size());
		
//		for (Point3d[] archetype: archetypes) {
//			System.out.println(Arrays.toString(archetype));
//		}
	
	    sc.stop();
	    sc.close();
	}
	
}
