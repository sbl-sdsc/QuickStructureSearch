package org.rcsb.project5;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

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
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.clustering.PowerIterationClustering;
import org.apache.spark.mllib.clustering.PowerIterationClusteringModel;
import org.biojava.nbio.structure.symmetry.geometry.SuperPositionQCP;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.hadoop.io.SimplePolymerType;
import org.rcsb.utils.BlastClustReader;

import scala.Tuple2;
import scala.Tuple3;

/**
 * This class clusters protein chains in 100% sequence identity clusters by structural similarity.
 * 
 * The input to the application consists of ten command line arguments:
 * 
 * <Hadoop sequence file name> <Output file name> <start index of first cluster> <end index for last cluster> <Maximum RMSD for structural clusters> <Interval added onto the maximum RMSD> <Number of intervals added> <Gap penalty> <Hole penalty> <Threshold output File Name>
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
		int option = Integer.parseInt(args[0]); // option 0 = print clusters; option 1 = threshold analysis; 
		// option 2 = print all clusters; option 3 = all threshold analysis
		String hadoopSequenceFileName = args[1];
		String outputFileName = args[2];
		int startCluster = Integer.parseInt(args[3]);
		int endCluster = Integer.parseInt(args[4]);
		double maxRmsd = Double.parseDouble(args[5]);
		double interval = Double.parseDouble(args[6]);
		int numInterval = Integer.parseInt(args[7]);
		double gapPenalty = Double.parseDouble(args[8]);
		double holePenalty = Double.parseDouble(args[9]);

		SequenceToStructureClusterer clusterer = new SequenceToStructureClusterer();
		if(option == 0) 
			clusterer.printClusters(hadoopSequenceFileName, startCluster, endCluster, outputFileName, maxRmsd, gapPenalty, holePenalty);
		else if(option == 1)
			clusterer.threshold(hadoopSequenceFileName, startCluster, endCluster, outputFileName, maxRmsd, interval, numInterval);
		else if(option == 2)
			clusterer.printAllClusters(hadoopSequenceFileName, outputFileName, maxRmsd, gapPenalty, holePenalty);
		else if(option == 3)
			clusterer.allThreshold(hadoopSequenceFileName, outputFileName, maxRmsd, interval, numInterval);
		else
			clusterer.runPIC(hadoopSequenceFileName, startCluster, endCluster, outputFileName, 2, gapPenalty, holePenalty);
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
	private void printClusters(String hadoopSequenceFileName, int startCluster, int endCluster, String outputFileName, double maxRmsd, double gapPenalty, double holePenalty) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();

		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc)
				.filter(c -> (c._1 >= startCluster && c._1 <= endCluster)).cache();
		//		System.out.println("Cluster pairs: " + clusterPairs.count());

		// create a list of chain ids needed for calculation
		List<String> chainIds = clusterPairs.values().collect();
		//		System.out.println("Chains: " + chainIds.size());

		// read <PdbId.chainId, SimplePolymerChain> pairs;
		JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
				//		.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(t -> (chainIds.contains(t._1)));
		//		System.out.println("chain count: " + chains.count());
		//		System.exit(-1);

		// broadcast <PdbId.chainId, SimplePolymerChain> pairs
		final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
		chains.unpersist();

		// group by cluster id
		JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();
		//		System.out.println("Clusters: " + clusters.count());

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
				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), list.get(0)._2[1].intValue(), SPCList, null, gapPenalty, holePenalty));
			} else {
				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList, null, gapPenalty, holePenalty));
			}
		}

		for(Cluster c: strClusterList) {
			c.findRepChain();
		}

		/*		for(Cluster c: strClusterList) {
			System.out.println(c);
		}*/

		for(Cluster c: strClusterList) {
			writeToCsv(writer, c);
		}

		writer.close();
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	private void threshold(String hadoopSequenceFileName, int startCluster, int endCluster, String outputFileName, double maxRmsd, double interval, int numInterval) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();

		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc)
				.filter(c -> (c._1 >= startCluster && c._1 <= endCluster)).cache();
		//		System.out.println("Cluster pairs: " + clusterPairs.count());

		// create a list of chain ids needed for calculation
		List<String> chainIds = clusterPairs.values().collect();
		//		System.out.println("Chains: " + chainIds.size());

		// read <PdbId.chainId, SimplePolymerChain> pairs;
		JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
				//		.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(t -> (chainIds.contains(t._1)));
		//		System.out.println("chain count: " + chains.count());
		//		System.exit(-1);

		// broadcast <PdbId.chainId, SimplePolymerChain> pairs
		final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
		chains.unpersist();

		// group by cluster id
		JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();
		//		System.out.println("Clusters: " + clusters.count());



		/*		for(Cluster c: strClusterList) {
			System.out.println(c);
		}*/



		/*for(List<Tuple2<String, Integer[]>> list: structuralClusterList) {
			writeToCsv(writer, list);
		}*/
		PrintWriter writerTest = new PrintWriter(outputFileName);
		writerTest.println("Threshold,Structural Clusters,Sequence Clusters");

		for (int i = 0; i <= numInterval; i++) {
			// loop through sequence clusters
			double threshold = maxRmsd + i*interval;

			writerTest.print(threshold);
			writerTest.print(",");

			List<List<Tuple2<String, Integer[]>>> structuralClusterList = clusters
					.map(t -> new StructureClusterer(chainMap, threshold)
					.getStructuralClusters(t)).collect();

			structuralClusterList = splitStrCluster(structuralClusterList);
			List<Cluster> strClusterList = new ArrayList<Cluster>();
			Map<String, SimplePolymerChain> map = chainMap.getValue();
			for (List<Tuple2<String, Integer[]>> list : structuralClusterList) {
				List<Tuple2<String, SimplePolymerChain>> SPCList = new ArrayList<Tuple2<String, SimplePolymerChain>>();
				for (Tuple2<String, Integer[]> tuple : list) {
					SPCList.add(new Tuple2<String, SimplePolymerChain>(tuple._1, map.get(tuple._1)));
				}
				if (list.get(0)._2[1] != null) {
					strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), list.get(0)._2[1].intValue(), SPCList));
				} else {
					strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList));
				}
			}
			int sequenceClusters = 0;
			int structuralClusters = 0;
			for (Cluster c : strClusterList) {
				if(c.getStrClusterId() >= 1) {
					structuralClusters++;
					if(c.getStrClusterId() == 1) {
						sequenceClusters++;
					}
				}
			}
			writerTest.print(structuralClusters);
			writerTest.print(",");
			writerTest.print(sequenceClusters);
			writerTest.println();
		}
		writerTest.flush();
		writerTest.close();
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	private void printAllClusters(String hadoopSequenceFileName, String outputFileName, double maxRmsd, double gapPenalty, double holePenalty) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();
		int numSeqClusters = getNumSeqClusters();
		PrintWriter writer = new PrintWriter(outputFileName);
		writer.println("SequenceClusterNumber, StructureClusterNumber, StructuralClusterSize, RepresentativeChain, OtherChains");

		int count = 0;

		for(int n = 0; n < numSeqClusters; n++) {
			int startCluster = n;
			int endCluster = n + count;
			//read <sequence cluster id, PdbId.chainId> pairs
			JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc)
					.filter(c -> ((c._1 >= startCluster) && (c._1 <= endCluster)))
					.cache();
			//		System.out.println("Cluster pairs: " + clusterPairs.count());

			// create a list of chain ids needed for calculation
			List<String> chainIds = clusterPairs.values().collect();
			//		System.out.println("Chains: " + chainIds.size());

			// read <PdbId.chainId, SimplePolymerChain> pairs;
			JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
					//		.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
					.filter(t -> (chainIds.contains(t._1)));
			//		System.out.println("chain count: " + chains.count());
			//		System.exit(-1);

			// broadcast <PdbId.chainId, SimplePolymerChain> pairs
			final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
			chains.unpersist();

			// group by cluster id
			JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();
			//		System.out.println("Clusters: " + clusters.count());

			//writer.println("PdbId.ChainId, SequenceClusterNumber, StructureClusterNumber");
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
					strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), list.get(0)._2[1].intValue(), SPCList, null, gapPenalty, holePenalty));
				} else {
					strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList, null, gapPenalty, holePenalty));
				}
			}

			for(Cluster c: strClusterList) {
				c.findRepChain();
			}

			/*		for(Cluster c: strClusterList) {
			System.out.println(c);
		}*/

			for(Cluster c: strClusterList) {
				writeToCsv(writer, c);
			}
			n += count;
			count ++;
		}
		writer.close();
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	private void allThreshold(String hadoopSequenceFileName, String outputFileName, double maxRmsd, double interval, int numInterval) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();
		int numSeqClusters = getNumSeqClusters();
		PrintWriter writerTest = new PrintWriter(outputFileName);

		int count = 0;

		writerTest.println("Threshold,Structural Clusters,Sequence Clusters");

		List<List<Double>> thresholds = new ArrayList<List<Double>> ();
		for(int x = 0; x <= numInterval; x ++) {
			List<Double> tempList = new ArrayList<Double> ();
			tempList.add((double)0);
			tempList.add((double)0);
			thresholds.add(tempList);
		}

		for(int n = 0; n < numSeqClusters; n++) {
			int startCluster = n;
			int endCluster = n + count;

			//read <sequence cluster id, PdbId.chainId> pairs
			JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc)
					.filter(c -> ((c._1 >= startCluster) && (c._1 <= endCluster)))
					.cache();
			//		System.out.println("Cluster pairs: " + clusterPairs.count());

			// create a list of chain ids needed for calculation
			List<String> chainIds = clusterPairs.values().collect();
			//		System.out.println("Chains: " + chainIds.size());

			// read <PdbId.chainId, SimplePolymerChain> pairs;
			JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
					//		.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
					.filter(t -> (chainIds.contains(t._1)));
			//		System.out.println("chain count: " + chains.count());
			//		System.exit(-1);

			// broadcast <PdbId.chainId, SimplePolymerChain> pairs
			final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
			chains.unpersist();

			// group by cluster id
			JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();
			//		System.out.println("Clusters: " + clusters.count());



			/*		for(Cluster c: strClusterList) {
			System.out.println(c);
		}*/



			/*for(List<Tuple2<String, Integer[]>> list: structuralClusterList) {
			writeToCsv(writer, list);
		}*/

			for (int i = 0; i <= numInterval; i++) {
				// loop through sequence clusters
				double threshold = maxRmsd + i*interval;

				int structuralClusters = 0;
				int sequenceClusters = 0;

				List<List<Tuple2<String, Integer[]>>> structuralClusterList = clusters
						.map(t -> new StructureClusterer(chainMap, threshold)
						.getStructuralClusters(t)).collect();

				structuralClusterList = splitStrCluster(structuralClusterList);
				List<Cluster> strClusterList = new ArrayList<Cluster>();
				Map<String, SimplePolymerChain> map = chainMap.getValue();
				for (List<Tuple2<String, Integer[]>> list : structuralClusterList) {
					List<Tuple2<String, SimplePolymerChain>> SPCList = new ArrayList<Tuple2<String, SimplePolymerChain>>();
					for (Tuple2<String, Integer[]> tuple : list) {
						SPCList.add(new Tuple2<String, SimplePolymerChain>(tuple._1, map.get(tuple._1)));
					}
					if (list.get(0)._2[1] != null) {
						strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), list.get(0)._2[1].intValue(), SPCList));
					} else {
						strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList));
					}
				}
				for (Cluster c : strClusterList) {
					if(c.getStrClusterId() >= 1) {
						structuralClusters++;
						if(c.getStrClusterId() == 1) {
							sequenceClusters++;
						}
					}
				}
				// System.out.println(thresholds.get(i).get(0));
				// System.out.println(thresholds.get(i).get(1));
				thresholds.get(i).set(0, (thresholds.get(i).get(0) + structuralClusters));
				thresholds.get(i).set(1, (thresholds.get(i).get(1) + sequenceClusters));
			}
			n += count;
			count ++;
		}
		for(int num = 0; num <= numInterval; num ++) {
			double threshold = maxRmsd + num*interval;
			writerTest.print(threshold);
			writerTest.print(",");
			writerTest.print(thresholds.get(num).get(0));
			writerTest.print(",");
			writerTest.print(thresholds.get(num).get(1));
			writerTest.println();
		}
		writerTest.flush();
		writerTest.close();
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
	
	private void runPIC(String hadoopSequenceFileName, int startCluster, int endCluster, String outputFileName, int groupNumber, double gapPenalty, double holePenalty) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();

		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc)
				.filter(c -> (c._1 >= startCluster && c._1 <= endCluster)).cache();
		//		System.out.println("Cluster pairs: " + clusterPairs.count());

		// create a list of chain ids needed for calculation
		List<String> chainIds = clusterPairs.values().collect();
		//		System.out.println("Chains: " + chainIds.size());

		// read <PdbId.chainId, SimplePolymerChain> pairs;
		JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
				//		.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(t -> (chainIds.contains(t._1)));
//		System.out.println("chain count: " + chains.count());

		// broadcast <PdbId.chainId, SimplePolymerChain> pairs
		final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
		chains.unpersist();

		// group by cluster id
		JavaPairRDD<Integer, Iterable<String>> clusters = clusterPairs.groupByKey();

		PrintWriter writer = new PrintWriter(outputFileName);
		//writer.println("PdbId.ChainId, SequenceClusterNumber, StructureClusterNumber");
		writer.println("SequenceClusterNumber, StructureClusterNumber, StructuralClusterSize, RepresentativeChain, OtherChains");

		List<List<Tuple2<String, Integer[]>>> structuralClusterList = new ArrayList<List<Tuple2<String, Integer[]>>>();
		List<Tuple2<Integer, List<Tuple2<String, SimplePolymerChain>>>> sequenceClusterList = clusters.map(t -> new StructureClusterer(chainMap).getSequenceCluster(t)).collect();
		for(Tuple2<Integer, List<Tuple2<String, SimplePolymerChain>>> sequenceCluster: sequenceClusterList) {
			int group = groupNumber;
			
			List<Tuple2<String, Integer[]>> nullCluster = new ArrayList<Tuple2<String, Integer[]>>();
			for(int i = sequenceCluster._2.size() - 1; i >= 0; i--) {
				if(sequenceCluster._2.get(i)._2 == null) {
					Tuple2<String, SimplePolymerChain> nullSPC = sequenceCluster._2.remove(i);
					nullCluster.add(new Tuple2<String, Integer[]>(nullSPC._1, new Integer[]{sequenceCluster._1, 0}));
				}
			}
			if(!nullCluster.isEmpty()) {
				structuralClusterList.add(nullCluster);
			}
			if(group > sequenceCluster._2.size()) {
				group = sequenceCluster._2.size();
			}
			if(group != 0) {
				List<Tuple3<Long, Long, Double>> affMatrix = getAffinityMatrix(sequenceCluster._2);
				JavaRDD<Tuple3<Long, Long, Double>> similarities = sc.parallelize(affMatrix);

				PowerIterationClustering pic = new PowerIterationClustering()
				.setK(group)
				.setMaxIterations(10);
				PowerIterationClusteringModel model = pic.run(similarities);

				List<List<Tuple2<String, Integer[]>>> newStructuralClusters = new ArrayList<List<Tuple2<String, Integer[]>>>();
				for(int i = 0; i < group; i++) {
					newStructuralClusters.add(new ArrayList<Tuple2<String, Integer[]>>());
				}
				
				for (PowerIterationClustering.Assignment a: model.assignments().toJavaRDD().collect()) {
					newStructuralClusters.get(a.cluster()).add(new Tuple2<String, Integer[]>(sequenceCluster._2.get((int) a.id())._1, new Integer[]{sequenceCluster._1, a.cluster() + 1}));
					//System.out.println(a.id() + " -> " + a.cluster());
				}
				for(List<Tuple2<String, Integer[]>> cluster: newStructuralClusters) {
					structuralClusterList.add(cluster);
				}
				//structuralClusterList.addAll(newStructuralClusters);
			}
		}

		// loop through sequence clusters
		//		List<List<Tuple2<String, Integer[]>>> structuralClusterList = clusters.map(t -> new StructureClusterer(chainMap, maxRmsd).getStructuralClusters(t)).collect();
//		structuralClusterList = splitStrCluster(structuralClusterList);
		List<Cluster> strClusterList = new ArrayList<Cluster>();
		Map<String, SimplePolymerChain> map = chainMap.getValue();
		for(List<Tuple2<String, Integer[]>> list: structuralClusterList) {
			List<Tuple2<String, SimplePolymerChain>> SPCList = new ArrayList<Tuple2<String, SimplePolymerChain>>();
			for(Tuple2<String, Integer[]> tuple: list) {
				SPCList.add(new Tuple2<String, SimplePolymerChain>(tuple._1, map.get(tuple._1)));
			}

			if (list.isEmpty()) {
			}else if (list.get(0)._2[1] != null) {
				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), list.get(0)._2[1].intValue(), SPCList, null, gapPenalty, holePenalty));
			} else {
				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList, null, gapPenalty, holePenalty));
			}
		}

		for(Cluster c: strClusterList) {
			c.findRepChain();
		}

		/*		for(Cluster c: strClusterList) {
			System.out.println(c);
		}*/

		for(Cluster c: strClusterList) {
			writeToCsv(writer, c);
		}

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
				//		.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(t -> t._2 != null)
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

	private static int getNumSeqClusters() {
		int sequenceIdentity = 100;
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);

		Map<String, Integer> clusterMap = reader.getPdbChainIdClusterMap();
		return clusterMap.size();
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
	
	public List<Tuple3<Long, Long, Double>> getAffinityMatrix(List<Tuple2<String, SimplePolymerChain>> seqCluster) {
		LongestCommonSubstring lcs = new LongestCommonSubstring();

		List<Tuple3<Long, Long, Double>> affMatrix = new ArrayList<Tuple3<Long, Long, Double>>();
		List<Integer> startEnd = null;

		for(int outer = 0; outer < seqCluster.size() - 1; outer++) {
			for(int inner = 1; inner < seqCluster.size(); inner++) {
				startEnd = lcs.longestCommonSubstring(
						seqCluster.get(outer)._2.getSequence(), seqCluster.get(inner)._2.getSequence());
				affMatrix.add(new Tuple3<Long, Long, Double>(new Long(outer), new Long(inner), (Double)getcRmsd(seqCluster.get(outer), seqCluster.get(inner), startEnd.get(0), startEnd.get(1), startEnd.get(2), startEnd.get(3))));
			}
		}
		return affMatrix;
	}

	/**
	 * Returns the cRMSD value of a specific portion of the Point3D arrays of the two chains given the starting and ending positions of the specific portions
	 * @param chain1
	 * @param chain2
	 * @param start1
	 * @param end1
	 * @param start2
	 * @param end2
	 * @return qcp.getRmsd()
	 */
	public double getcRmsd(Tuple2<String, SimplePolymerChain> chain1, Tuple2<String, SimplePolymerChain> chain2, int start1, int end1, int start2, int end2) {
		SuperPositionQCP qcp = new SuperPositionQCP();
		int s1 = start1;
		int e1 = end1;
		int s2 = start2;
		List<Point3d> coordinates1 = new ArrayList<Point3d>();
		List<Point3d> coordinates2 = new ArrayList<Point3d>();
		for(int n = 0; n < (e1 - s1); n++) {
			Point3d temp1 = chain1._2.getCoordinates()[s1 + n];
			Point3d temp2 = chain2._2.getCoordinates()[s2 + n];
			if((temp1 != null)&&(temp2 != null)) {
				coordinates1.add(temp1);
				coordinates2.add(temp2);
			}
		}
		Point3d[] seq1 = new Point3d[coordinates1.size()];
		Point3d[] seq2 = new Point3d[coordinates2.size()];
		for(int count = 0; count < seq1.length; count ++) {
			seq1[count] = coordinates1.get(count);
			seq2[count] = coordinates2.get(count);
		}
		qcp.set(seq1, seq2);
		return qcp.getRmsd();
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
