package org.rcsb.project5;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.hadoop.io.SimplePolymerType;
import org.rcsb.utils.BlastClustReader;

import scala.Tuple2;

/**
 * This class clusters protein chains in 100% sequence identity clusters by structural similarity.
 * 
 * The input to the application consists of three command line arguments:
 * 
 * <Hadoop sequence file name> <start index of first cluster> <end index for last cluster>
 * 
 * To run a small example, use a small range of clusters, i.e. 5000 - 5010.
 * 
 * @author  Peter Rose
 */
public class SequenceToStructureClusterer {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String hadoopSequenceFileName = args[0];
		int startCluster = Integer.parseInt(args[1]);
        int endCluster = Integer.parseInt(args[2]);

        SequenceToStructureClusterer clusterer = new SequenceToStructureClusterer();
        clusterer.run(hadoopSequenceFileName, startCluster, endCluster);
	}

	/**
	 * Performs structural clustering for the sequence clusters within a specified cluster index range.
	 * 
	 * @param hadoopSequenceFileName
	 * @param startCluster index of first cluster
	 * @param endCluster index of last cluster
	 */
	private void run(String hadoopSequenceFileName, int startCluster, int endCluster) {
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();
		
		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc)
				.filter(c -> (c._1 >= startCluster && c._1 <= endCluster)).cache();
		System.out.println("Cluster pairs: " + clusterPairs.count());
		
		// create a list of chain ids needed for calculation
	    List<String> chainIds = clusterPairs.map(c -> c._2).collect();
		System.out.println("Chains: " + chainIds.size());
		
		// read <PdbId.chainId, SimplePolymerChain> pairs;
		JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
				.filter(t -> (chainIds.contains(t._1)));
		System.out.println("chain count: " + chains.count());
		
	    // broadcast <PdbId.chainId, SimplePolymerChain> pairs
		final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
		chains.unpersist();
				
		// group by cluster id
		JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();
		System.out.println("Clusters: " + clusters.count());
		
		// loop through sequence clusters
		clusters.foreach(new StructureClusterer(chainMap));
		
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	/**
	 * Gets pairs of chain id, simple polymer chains for proteins from a Hadoop sequence file. 
	 * 
	 * @param hadoopSequenceFileName file containing simple polymer chains
	 * 
	 * @param sc
	 * @return pairs of chain id, simple polymer chain
	 */
	private static JavaPairRDD<String, SimplePolymerChain> getPolymerChains(String hadoopSequenceFileName, JavaSparkContext sc) {
		JavaPairRDD<String, SimplePolymerChain> chains = sc
						.sequenceFile(hadoopSequenceFileName, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
						.mapToPair(new HadoopToSimpleChainMapper())
						.filter(t -> t._2.isProtein());
		return chains;
	}
	
	/**
	 * Gets pairs of sequence cluster identifier, chain id.
	 * 
	 * @param sc
	 * @return pairs of sequence cluster identifier, chain id
	 */
	private static JavaPairRDD<Integer, String> getSequenceClusters(JavaSparkContext sc) {
		int sequenceIdentity = 100;
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);

		Map<String, Integer> clusterMap = reader.getPdbChainIdClusterMap();
		List<Tuple2<Integer, String>> list = new ArrayList<>(clusterMap.size());
		for (Entry<String, Integer> entity: clusterMap.entrySet()) {
			list.add(new Tuple2<Integer,String>(entity.getValue(), entity.getKey()));
		}
		JavaPairRDD<Integer, String> clusterPairs = sc.parallelizePairs(list, NUM_THREADS*NUM_TASKS_PER_THREAD).cache();
		
		return clusterPairs;
	}
	
	/**
	 * Collects a JavaPairRDD<K, V> to a Map<K,V>
	 * @param pairRdd the JavaPairRDD<K, V> to convert
	 * @return Map<K, V>
	 */
	private static <K, V> Map<K, V> collectAsMap(JavaPairRDD<K, V> pairRdd) {
		List<Tuple2<K, V>> tuples = pairRdd.collect();
		Map<K, V> map = new HashMap<>(tuples.size());
		
		for (Tuple2<K, V> t: tuples) {
			map.put(t._1,  t._2);
		}
		
		return map;
	}
	
	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(SequenceToStructureClusterer.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
//				.set("spark.kryoserializer.buffer.max", "1000m")
//				.set("spark.driver.maxResultSize", "1000m");
		
		conf.registerKryoClasses(new Class[]{SimplePolymerChain.class, SimplePolymerType.class});
		
		return new JavaSparkContext(conf);
	}
}
