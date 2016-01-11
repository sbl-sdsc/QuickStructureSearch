package org.rcsb.hadoop.io;

import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.Comparator;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.compress.IntegerDeltaZigzagVariableByte;
import org.rcsb.compress.IntegerToByteTransform;
import org.rcsb.structuralSimilarity.RadiusOfGyrationMapper;

import scala.Tuple2;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file with
 * PDBId.ChainId/SimplePolymerChain pairs. This demo program calculates the radius of gyration
 * for each polymer chain and calculates some basic statistics. It also prints out the polymer 
 * chains with the top 10% highest radius of gyration.
 * 
 * @author  Peter Rose
 */
public class HadoopSequenceFileRGyration implements Serializable {    
	private static final long serialVersionUID = 1L;
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
//	private static final Logger logger = LoggerFactory.getLogger(HadoopSequenceFileRGyration.class);
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0]; // Hadoop sequence input file name

		// This is the default 2 line structure for Spark applications
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(HadoopSequenceFileRGyration.class.getSimpleName())
				.set("spark.ui.showConsoleProgress", "false")
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		
		Class<?>[] classes = {SimplePolymerChain.class, IntegerToByteTransform.class, IntegerDeltaZigzagVariableByte.class, RadiusOfGyrationMapper.class};
		conf.registerKryoClasses(classes);
		
		JavaSparkContext sc = new JavaSparkContext(conf);

		long start = System.nanoTime();
		
    	// read sequence file and map sequence length to an RDD
		JavaPairRDD<String, Float> chains = sc
				.sequenceFile(path, Text.class, SimplePolymerChain.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1.toString(), t._2.getCoordinates()) )	
				.mapValues(new RadiusOfGyrationMapper())
				.cache();
		
		// map radius of gyration to a new JavaRDD to calculate various statistics
		JavaRDD<Float> radiusOfGyration = chains
				.map(t -> t._2)
				.cache();
		
		// calculate basis statistics on radius of gyration
		long count = chains.count();
		float min = radiusOfGyration.min(Comparator.naturalOrder());
		float max = radiusOfGyration.max(Comparator.naturalOrder());	
		float mean = radiusOfGyration.reduce((a,b) -> a + b)/count;
		
		// find chains with the top 10% of radius of gyration
		float threshold = max - (max - min)/10;
		chains
		        .filter(t -> t._2 > threshold)
		        .foreach(t -> System.out.println(t));
		
		long end = System.nanoTime();
		
		// Shutdown Spark
		sc.stop();
		sc.close();
				
        System.out.println("Total chains               : " + count);	
		System.out.println("Mean radius of gyration    : " + mean);
		System.out.println("Minimum radius of gyration : " + min);
		System.out.println("Maximum radius of gyration : " + max);
		
		System.out.println("Time: " + (end - start)/1E9 + " sec.");
	}
}
