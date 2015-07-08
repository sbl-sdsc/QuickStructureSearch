package org.rcsb.project5;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
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
 * The input to the application consists of five command line arguments:
 * 
 * <Hadoop sequence file name> <Output file name> <start index of first cluster> <end index for last cluster> <Maximum RMSD for structural clusters>
 * 
 * To run a small example, use a small range of clusters, i.e. 5000 - 5010.
 * 
 * @author Peter Rose
 * @author Justin Li
 * @author Joe Sun
 */
public class SequenceToStructureClusterer {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String hadoopSequenceFileName = args[0];
		String outputFileName = args[1];
		int startCluster = Integer.parseInt(args[2]);
		int endCluster = Integer.parseInt(args[3]);
		double maxRmsd = Double.parseDouble(args[4]);

		SequenceToStructureClusterer clusterer = new SequenceToStructureClusterer();
		clusterer.run(hadoopSequenceFileName, startCluster, endCluster, outputFileName, maxRmsd);
	}

	/**
	 * Performs structural clustering for the sequence clusters within a specified cluster index range.
	 * 
	 * @param hadoopSequenceFileName
	 * @param outputFileName
	 * @param startCluster index of first cluster
	 * @param endCluster index of last cluster
	 * @param maxRmsd
	 */
	private void run(String hadoopSequenceFileName, int startCluster, int endCluster, String outputFileName, double maxRmsd) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();

		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc)
				.filter(c -> (c._1 >= startCluster && c._1 <= endCluster)).cache();
		System.out.println("Cluster pairs: " + clusterPairs.count());

		// create a list of chain ids needed for calculation
		List<String> chainIds = clusterPairs.values().collect();
		System.out.println("Chains: " + chainIds.size());

		// read <PdbId.chainId, SimplePolymerChain> pairs;
		JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
				.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(t -> (chainIds.contains(t._1)));
		System.out.println("chain count: " + chains.count());

		// broadcast <PdbId.chainId, SimplePolymerChain> pairs
		final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
		chains.unpersist();

		// group by cluster id
		JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();
		System.out.println("Clusters: " + clusters.count());

		PrintWriter writer = new PrintWriter(outputFileName);
		//writer.println("PdbId.ChainId, SequenceClusterNumber, StructureClusterNumber");
		writer.println("SequenceClusterNumber, StructureClusterNumber, StructuralClusterSize, RepresentativeChain, OtherChains");
		// loop through sequence clusters
		List<List<Tuple2<String, Integer[]>>> structuralClusterList = clusters.map(t -> new StructureClusterer(chainMap, maxRmsd).getStructuralClusters(t)).collect();
		structuralClusterList = splitStrCluster(structuralClusterList);
		List<Cluster> strClusterList = new ArrayList<Cluster>();
		Map<String, SimplePolymerChain> map = chainMap.getValue();
		for(List<Tuple2<String, Integer[]>> list: structuralClusterList) {
			List<Tuple2<String, SimplePolymerChain>> SPCList = new ArrayList<Tuple2<String, SimplePolymerChain>>();
			for(Tuple2<String, Integer[]> tuple: list) {
				SPCList.add(new Tuple2<String, SimplePolymerChain>(tuple._1, map.get(tuple._1)));
			}

			if (list.get(0)._2[1] != null) {
				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), list.get(0)._2[1].intValue(), SPCList, null));
			} else {
				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList, null));
			}
		}
		
		for(Cluster c: strClusterList) {
			c.findRepChain();
		}

		for(Cluster c: strClusterList) {
			System.out.println(c);
		}

		for(Cluster c: strClusterList) {
			writeToCsv(writer, c);
		}

		/*for(List<Tuple2<String, Integer[]>> list: structuralClusterList) {
			writeToCsv(writer, list);
		}*/

		writer.close();
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
				.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
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
	
	private List<List<Tuple2<String, Integer[]>>> splitStrCluster(List<List<Tuple2<String, Integer[]>>> structuralClusterList) {
		List<List<Tuple2<String, Integer[]>>> strList = new ArrayList<List<Tuple2<String, Integer[]>>>();
		for (List<Tuple2<String, Integer[]>> tupleList: structuralClusterList) {
			for(int i = 0; i < tupleList.size(); i++) {
				boolean unique = true;
				for(List<Tuple2<String, Integer[]>> strCluster: strList) {
					//if the sequence and structural cluster numbers are the same
					if((strCluster.get(0)._2[0].intValue() == tupleList.get(i)._2[0].intValue())) {
						if(strCluster.get(0)._2[1] == null || tupleList.get(i)._2[1] == null) {
							if(strCluster.get(0)._2[1] == tupleList.get(i)._2[1]) {
								strCluster.add(tupleList.get(i));
								unique = false;
								break;
							}
						}
						else if (strCluster.get(0)._2[1].intValue() == tupleList.get(i)._2[1].intValue()) {
							strCluster.add(tupleList.get(i));
							unique = false;
							break;
						}
					}
				}
				if(unique) {
					List<Tuple2<String, Integer[]>> newStrCluster = new ArrayList<Tuple2<String, Integer[]>>();
					newStrCluster.add(tupleList.get(i));
					strList.add(newStrCluster);
				}
			}
		}
		return strList;
	}

	/**
	 * Writes chain id, sequence cluster, and structural cluster to a csv file
	 * @param writer
	 * @param list
	 */
	private static void writeToCsv(PrintWriter writer, List<Tuple2<String, Integer[]>> list) {
		for (Tuple2<String, Integer[]> t : list) {
			writer.print(t._1);
			for (Integer f: t._2) {
				writer.print(",");
				writer.print(f);
			}
			writer.println();
		}
		writer.flush();
	}

	/**
	 * Writes chain id, sequence cluster, and structural cluster to a csv file
	 * @param writer
	 * @param list
	 */
	private static void writeToCsv(PrintWriter writer, Cluster list) {
		writer.print(list.getSeqClusterId());
		writer.print(",");
		writer.print(list.getStrClusterId());
		writer.print(",");
		writer.print(list.size());
		writer.print(",");
		if (list.getRepChain() == null) {
			writer.print("null");
		} else {
			writer.print(list.getRepChain()._1);
		}
		for (Tuple2<String,SimplePolymerChain> t : list.getStrCluster()) {
			if (list.getRepChain() == null || !t._1.equals(list.getRepChain()._1)) {
				writer.print(",");
				writer.print(t._1);
			}
		}
		writer.println();
		writer.flush();
	}
}
