package org.rcsb.project5;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralAlignment.SuperPositionQCP;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.SeqToChainMapper;
import org.rcsb.utils.BlastClustReader;

import scala.Tuple2;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class FindClusters {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(FindClusters.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
				.set("spark.kryoserializer.buffer.max", "1000m");
		
		JavaSparkContext sc = new JavaSparkContext(conf);
		

		long start = System.nanoTime();
		
		// read sequence file and map to PdbId.chainId, SimplePolymerChain pairs
		JavaPairRDD<String, SimplePolymerChain> chains = sc
						.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
	//					.sample(false, 0.01, 123)
						.mapToPair(new HadoopToSimpleChainMapper()) // convert input to <pdbId.chainId, SimplePolymerChain> pairs
						.filter(t -> t._2.isProtein())
						.cache();

		System.out.println("number of chains: " + chains.count());
		int sequenceIdentity = 100;
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);
		
		List<List<String>> clusters = reader.getPdbChainIdClusters();
		
		for (List<String> cluster: clusters) {
			if (cluster.size() == 6) {
				chains.filter(t -> cluster.contains(t._1));
                List<Tuple2<String, SimplePolymerChain>> list = chains.collect();
                for (Tuple2<String, SimplePolymerChain> t: list) {
                	System.out.println(t._1 + ", " + t._2);
                }
                break;
			}
		}
		
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
}
