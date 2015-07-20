package org.rcsb.project5;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.hadoop.io.SimplePolymerType;
import org.rcsb.utils.BlastClustReader;

import scala.Tuple2;

public class SequenceClusterIdComparison {
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args) {
		
		int seqIdBig = 100;
		int seqIdSmall = 90;

		// initialize Spark		
		JavaSparkContext sc = getSparkContext();
		
		/*long startTest = System.nanoTime();
		List<Tuple2<Integer, List<Integer>>> list = getSeqIdClustering(90, sc);
		System.out.println(list);
		System.out.println("Time: " + (System.nanoTime() - startTest)/1E9 + " sec.");
		System.exit(-1);*/
		
		long start = System.nanoTime();

		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc, seqIdBig)
				.cache();

		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer,String> clusterPairs2 = getSequenceClusters(sc, seqIdSmall)
				.cache();

		// group by cluster id
		JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();

		// group by cluster id
		JavaPairRDD<Integer, Iterable<String>> clusters2 = clusterPairs2.groupByKey();

		int big = getNumSeqClusters(seqIdBig);
		int small = getNumSeqClusters(seqIdSmall);
		int score = big - small;
		
		System.out.println(big);
		System.out.println(small);
		System.out.println("The difference in cluster numbers is " + score);
//		System.exit(-1);
		
		List<Tuple2<Integer, Iterable<String>>> clustersList = clusters.map(t -> t).collect();
		
		List<Integer> dif = clusters2.map(t -> getGroupDifference(t._2, clustersList)).collect();
		
		int sum = 0;
		
		for(Integer i: dif) {
			System.out.println(i);
			sum += i;
		}
		
		System.out.println(big);
		System.out.println(small);
		System.out.println("The difference in cluster numbers is " + score);
		System.out.println("The value is " + sum);
		System.out.println("The difference between score and sum is " + (score - sum));
		
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
		
	}
	
	public static List<Tuple2<Integer, List<Integer>>> getSeqIdClustering(int seqId, JavaSparkContext sc) {
				//read <sequence cluster id, PdbId.chainId> pairs
				JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc, 100)
						.cache();

				//read <sequence cluster id, PdbId.chainId> pairs
				JavaPairRDD<Integer,String> clusterPairs2 = getSequenceClusters(sc, seqId)
						.cache();

				// group by cluster id
				JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();

				// group by cluster id
				JavaPairRDD<Integer, Iterable<String>> clusters2 = clusterPairs2.groupByKey();

				
				List<Tuple2<Integer, Iterable<String>>> clustersList = clusters.map(t -> t).collect();
				
				List<Tuple2<Integer, List<Integer>>> dif = clusters2.map(t -> new Tuple2<Integer, List<Integer>>(t._1, getGroupDifferenceArray(t._2, clustersList))).collect();
				
				return dif;
	}
	
	private static int getGroupDifference(Iterable<String> list, List<Tuple2<Integer, Iterable<String>>> clusters) {
		List<Tuple2<Integer, Iterable<String>>> clusterList = clusters;
		List<Integer> clusterIdList = new ArrayList<Integer>();
		for(Tuple2<Integer, Iterable<String>> tuple: clusterList) {
			for(String s: list) {
				if(contains(tuple, s)) {
					clusterIdList.add(tuple._1);
					break;
				}
			}
		}
		return clusterIdList.size() - 1;
	}
	
	private static List<Integer> getGroupDifferenceArray(Iterable<String> list, List<Tuple2<Integer, Iterable<String>>> clusters) {
		List<Tuple2<Integer, Iterable<String>>> clusterList = clusters;
		List<Integer> clusterIdList = new ArrayList<Integer>();
		for(Tuple2<Integer, Iterable<String>> tuple: clusterList) {
			for(String s: list) {
				if(contains(tuple, s)) {
					clusterIdList.add(tuple._1);
					break;
				}
			}
		}
		return clusterIdList;
	}
	
	private static boolean contains(Tuple2<Integer, Iterable<String>> tuple, String s) {
		Iterator<String> iter = tuple._2().iterator();
		while(iter.hasNext()) {
			if(iter.next().equals(s)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Gets pairs of sequence cluster identifier, chain id.
	 * 
	 * @param sc
	 * @return pairs of sequence cluster identifier, chain id
	 */
	private static JavaPairRDD<Integer, String> getSequenceClusters(JavaSparkContext sc, int seqId) {
		int sequenceIdentity = seqId;
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);

		Map<String, Integer> clusterMap = reader.getPdbChainIdClusterMap();
		List<Tuple2<Integer, String>> list = new ArrayList<>(clusterMap.size());
		for (Entry<String, Integer> entity: clusterMap.entrySet()) {
			list.add(new Tuple2<Integer,String>(entity.getValue(), entity.getKey()));
		}
		JavaPairRDD<Integer, String> clusterPairs = sc.parallelizePairs(list, NUM_THREADS*NUM_TASKS_PER_THREAD).cache();

		return clusterPairs;
	}

	private static int getNumSeqClusters(int seqId) {
		int sequenceIdentity = seqId;
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);

		List<List<String>> clusters = reader.getPdbChainIdClusters();
		return clusters.size();
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
