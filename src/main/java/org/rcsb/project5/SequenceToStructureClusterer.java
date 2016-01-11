package org.rcsb.project5;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;

import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.clustering.PowerIterationClustering;
import org.apache.spark.mllib.clustering.PowerIterationClusteringModel;
import org.biojava.nbio.structure.symmetry.geometry.SuperPositionQCP;
import org.rcsb.hadoop.io.HadoopToSimplePolymerChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.hadoop.io.SimplePolymerType;
import org.rcsb.utils.BlastClustReader;

import scala.Tuple2;
import scala.Tuple3;

/**
 * This class clusters protein chains in clusters of a chosen sequence identity by structural similarity.
 * 
 * The input to the application consists of twelve command line arguments:
 * 
 * <Option number> <Hadoop sequence file name> <Output file name> <start index of first cluster> <end index for last cluster> <Maximum RMSD for structural clusters> <Interval added onto the maximum RMSD> <Number of intervals added> <Gap penalty> <Hole penalty> <% sequence identity> <CSV input file name>
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
		// option 4 = PIC; option 5 = altered sequence identity preparation
		String hadoopSequenceFileName = args[1];
		String outputFileName = args[2];
		int startCluster = Integer.parseInt(args[3]);
		int endCluster = Integer.parseInt(args[4]);
		double maxRmsd = Double.parseDouble(args[5]);
		double interval = Double.parseDouble(args[6]);
		int numInterval = Integer.parseInt(args[7]);
		double gapPenalty = Double.parseDouble(args[8]);
		double holePenalty = Double.parseDouble(args[9]);
		int sequenceIdentity = Integer.parseInt(args[10]);
		String outputFileName2 = args[11];
		//		String csvFileName = args[11];

		SequenceToStructureClusterer clusterer = new SequenceToStructureClusterer();
		if(option == 0) 
			clusterer.printClusters(hadoopSequenceFileName, startCluster, endCluster, outputFileName, maxRmsd, gapPenalty, holePenalty);
		else if(option == 1)
			clusterer.threshold(hadoopSequenceFileName, startCluster, endCluster, outputFileName, maxRmsd, interval, numInterval);
		else if(option == 2)
			clusterer.printAllClusters(hadoopSequenceFileName, outputFileName, maxRmsd, gapPenalty, holePenalty);
		else if(option == 3)
			clusterer.allThreshold(hadoopSequenceFileName, outputFileName, maxRmsd, interval, numInterval);
		else if(option == 4)
			clusterer.runPIC(hadoopSequenceFileName, startCluster, endCluster, outputFileName, 2, gapPenalty, holePenalty);
		else
			clusterer.alteredSeqIdentity(hadoopSequenceFileName, outputFileName, startCluster, endCluster, maxRmsd, gapPenalty, holePenalty, sequenceIdentity, outputFileName2);
	}

	/**
	 * Performs structural clustering for the 100% sequence clusters within a specified cluster index range.
	 * 
	 * @param hadoopSequenceFileName the location of the hadoop sequence file
	 * @param outputFileName the desired location of the csv file that will be created
	 * @param startCluster index of first cluster
	 * @param endCluster index of last cluster
	 * @param maxRmsd the maximum RMSD value between any two chains of the same structural cluster (default: 1.0)
	 * @param gapPenalty the penalty given to each gap when calculating for the representative chain (default: 1.0)
	 * @param holePenalty the penalty given to each hole when calculating for the representative chain (default: 0.1)
	 */
	private void printClusters(String hadoopSequenceFileName, int startCluster, int endCluster, String outputFileName, double maxRmsd, double gapPenalty, double holePenalty) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();

		JavaPairRDD<Integer, String> clusterPairs;
		//read <sequence cluster id, PdbId.chainId> pairs
		clusterPairs = getSequenceClusters(sc)
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
			} 
			//			else {
			//				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList, null, gapPenalty, holePenalty));
			//			}
		}

		for(Cluster c: strClusterList) {
			c.findRepChain();
		}

		quickSort(strClusterList, 0, strClusterList.size() - 1);

		//		int prevSeqId = 0;
		//		int seqClusterId;
		//		int startIndex = 0;
		//		for(int n = 0; n < strClusterList.size(); n ++) {
		//			seqClusterId = strClusterList.get(n).getSeqClusterId();
		//			if(seqClusterId != prevSeqId) {
		//				if(n - 1 - startIndex > 1) {
		//					quickSortByStrId(strClusterList, startIndex, n-1);
		//				}
		//				startIndex = n;
		//				prevSeqId = strClusterList.get(n).getSeqClusterId();
		//			}
		//			else if(n == strClusterList.size()-1) {
		//				if(n - startIndex > 1)
		//					quickSortByStrId(strClusterList, startIndex, n);
		//			}
		//		}

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
	 * Analyzes the effect that the maximum RMSD threshold has on the number of structural clusters created for a given range of 100% sequence clusters
	 * Note: the range of the thresholds which will be analyzed will range from maxRmsd to maxRmsd + (interval)*(numInterval), inclusive
	 * @param hadoopSequenceFileName the location of the hadoop sequence file
	 * @param startCluster the index of the first sequence cluster which will be analyzed
	 * @param endCluster the index of the last sequence cluster which will be analyzed
	 * @param outputFileName the desired location of the csv file that will be created
	 * @param maxRmsd the lowest maximum RMSD threshold value which will be analyzed
	 * @param interval the length of each interval
	 * @param numInterval the number of intervals
	 * @throws FileNotFoundException
	 */
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

	/**
	 * Performs structural clustering for all 100% sequence clusters
	 * @param hadoopSequenceFileName the location of the hadoop sequence file
	 * @param outputFileName the desired location of the csv file that will be created
	 * @param maxRmsd the maximum RMSD value between any two chains of the same structural cluster (default: 1.0)
	 * @param gapPenalty the penalty given to each gap when calculating for the representative chain (default: 1.0)
	 * @param holePenalty the penalty given to each hole when calculating for the representative chain (default: 0.1)
	 * @throws FileNotFoundException
	 */
	private void printAllClusters(String hadoopSequenceFileName, String outputFileName, double maxRmsd, double gapPenalty, double holePenalty) throws FileNotFoundException{
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();
		//	int numSeqClusters = getNumSeqClusters();
		PrintWriter writer = new PrintWriter(outputFileName);
		writer.println("SequenceClusterNumber, StructureClusterNumber, StructuralClusterSize, RepresentativeChain, OtherChains");

		//	int ends[] = {-1, 0, 100, 1000, 5000, 10000, 20000, 30000, 45000, numSeqClusters};
		int ends[] = {999, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};

		for(int count = 0; count < ends.length - 1; count ++) {
			int startCluster = ends[count] + 1;
			int endCluster = ends[count + 1];
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
				} 
				//				else {
				//					strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList, null, gapPenalty, holePenalty));
				//				}
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
		}
		writer.close();
		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	/**
	 * Analyzes the effect that the maximum RMSD threshold has on the number of structural clusters created all 100% sequence clusters
	 * Note: the range of the thresholds which will be analyzed will range from maxRmsd to maxRmsd + (interval)*(numInterval), inclusive
	 * @param hadoopSequenceFileName the location of the hadoop sequence file
	 * @param outputFileName the desired location of the csv file that will be created
	 * @param maxRmsd the lowest maximum RMSD threshold value which will be analyzed
	 * @param interval the length of each interval
	 * @param numInterval the number of intervals
	 * @throws FileNotFoundException
	 */
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
	 * Prints sequence clusters of a given sequence identity of a certain range, further broken 
	 * down into structural clusters, into a csv file
	 * 
	 * @param hadoopSequenceFileName the location of the hadoop sequence file
	 * @param outputFileName the desired location of the csv file that will be created
	 * @param startCluster the index of the first sequence cluster to be printed
	 * @param endCluster the index of the last sequence cluster to be printed
	 * @param maxRmsd the maximum RMSD value between chains belonging to the same structural cluster (default: 1.0)
	 * @param gapPenalty the penalty given to each gap when calculating for the representative chain (default: 1.0)
	 * @param holePenalty the penalty given to each hole when calculating for the representative chain (default: 0.1)
	 * @param sequenceIdentity the percent sequence similarity between members of the same sequence cluster
	 * @throws FileNotFoundException 
	 */
	private void alteredSeqIdentity(String hadoopSequenceFileName, String outputFileName, int startCluster, int endCluster, double maxRmsd, double gapPenalty, double holePenalty, int sequenceIdentity, String outputFileName2) throws FileNotFoundException {
		// initialize Spark		
		JavaSparkContext sc = getSparkContext();

		long start = System.nanoTime();
		List<Integer> seq100Clusters = SequenceClusterIdComparison.getSeq100Clusters(sequenceIdentity, sc, startCluster, endCluster);
		for(int n = 0; n < seq100Clusters.size(); n ++) {
			System.out.println(seq100Clusters.get(n));
		}
		//read <sequence cluster id, PdbId.chainId> pairs
		JavaPairRDD<Integer, String> clusterPairs = getSequenceClusters(sc, 100)
				.filter(c -> seq100Clusters.contains(c._1))
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
			} 
			//			else {
			//				strClusterList.add(new Cluster(list.get(0)._2[0].intValue(), 0, SPCList, null, gapPenalty, holePenalty));
			//			}
		}

		for(Cluster c: strClusterList) {
			c.findRepChain();
		}

		quickSort(strClusterList, 0, strClusterList.size() - 1);
		
		PrintWriter writer100 = new PrintWriter(outputFileName2);
		writer100.println("SequenceClusterNumber, StructureClusterNumber, StructuralClusterSize, RepresentativeChain, OtherChains");
		for(Cluster c: strClusterList) {
			writeToCsv(writer100, c);
		}

		writer100.close();

		List<List<Cluster>> clusterList = new ArrayList<List<Cluster>>();

		int prevSeqId = -1;
		int seqClusterId;
		int count = -1;

		for(int a = 0; a < strClusterList.size(); a ++) {
			seqClusterId = strClusterList.get(a).getSeqClusterId();
			if(seqClusterId != prevSeqId) {
				clusterList.add(new ArrayList<Cluster>());
				count ++;
				clusterList.get(count).add(strClusterList.get(a));
				prevSeqId = seqClusterId;
			}
			else {
				clusterList.get(count).add(strClusterList.get(a));
			}
		}

		for(int x = 0; x < clusterList.size(); x ++) {
			quickSort(clusterList.get(x), 0, clusterList.get(x).size() - 1);
		}

		//		BufferedReader reader = null;
		//		List<Integer> seq100Clusters = SequenceClusterIdComparison.getSeq100Clusters(sequenceIdentity, sc, startCluster, endCluster);
		//
		//		//read <sequence cluster id, PdbId.chainId> pairs
		//		JavaPairRDD<Integer,String> clusterPairs = getSequenceClusters(sc, 100)
		//				.filter(c -> seq100Clusters.contains(c._1))
		//				.cache();
		//		//		System.out.println("Cluster pairs: " + clusterPairs.count());
		//
		//		// create a list of chain ids needed for calculation
		//		List<String> chainIds = clusterPairs.values().collect();
		//		//		System.out.println("Chains: " + chainIds.size());
		//
		//		// read <PdbId.chainId, SimplePolymerChain> pairs;
		//		JavaPairRDD<String, SimplePolymerChain> chains = getPolymerChains(hadoopSequenceFileName, sc)
		//				//		.filter(new GapFilterSPC(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
		//				.filter(t -> (chainIds.contains(t._1)));
		//		//		System.out.println("chain count: " + chains.count());
		//		//		System.exit(-1);
		//
		//		// broadcast <PdbId.chainId, SimplePolymerChain> pairs
		//		final Broadcast<Map<String, SimplePolymerChain>> chainMap = sc.broadcast(collectAsMap(chains));
		//		chains.unpersist();
		//
		//		Map<String, SimplePolymerChain> map = chainMap.getValue();
		//
		//		List<List<Cluster>> clusterList = new ArrayList<List<Cluster>> (); 
		//
		//		int prevSeqClusterId = -1;
		//		int seqClusterId;
		//		int strClusterId;
		//		int clusterSize;
		//		int count = -1;
		//
		//		try {
		//			reader = new BufferedReader(new FileReader(csvFileName));
		//			String line = "";
		//
		//			while((line = reader.readLine()) != null) {
		//				String[] tokens = line.split(",");
		//				seqClusterId = Integer.parseInt(tokens[0]);
		//				if(seqClusterId != prevSeqClusterId) {
		//					count ++;
		//					clusterList.add(new ArrayList<Cluster> ());
		//					strClusterId = Integer.parseInt(tokens[1]);
		//					clusterSize = Integer.parseInt(tokens[2]);
		//					List<Tuple2<String, SimplePolymerChain>> strCluster = new ArrayList<Tuple2<String, SimplePolymerChain>> ();
		//					for(int n = 0; n < clusterSize; n ++) {
		//						String chainId = tokens[n + 3];
		//						strCluster.add(new Tuple2<String, SimplePolymerChain> (chainId, map.get(chainId)));
		//					}
		//					clusterList.get(count).add(new Cluster(seqClusterId, strClusterId, strCluster, null, gapPenalty, holePenalty));
		//					prevSeqClusterId = seqClusterId;
		//				}
		//				else {
		//					strClusterId = Integer.parseInt(tokens[1]);
		//					clusterSize = Integer.parseInt(tokens[2]);
		//					List<Tuple2<String, SimplePolymerChain>> strCluster = new ArrayList<Tuple2<String, SimplePolymerChain>> ();
		//					for(int n = 0; n < clusterSize; n ++) {
		//						String chainId = tokens[n + 3];
		//						strCluster.add(new Tuple2<String, SimplePolymerChain> (chainId, map.get(chainId)));
		//					}
		//					clusterList.get(count).add(new Cluster(seqClusterId, strClusterId, strCluster, null, gapPenalty, holePenalty));
		//				}
		//			}
		//		}
		//		catch (Exception e) {
		//			e.printStackTrace();
		//		}
		//		finally {
		//			try {
		//				reader.close();
		//			}
		//			catch (IOException e) {
		//				e.printStackTrace();
		//			}
		//		}
		//		for(int a = 0; a < clusterList.size(); a ++) {
		//			for(int b = 0; b < clusterList.get(a).size(); b ++) {
		//				System.out.print(clusterList.get(a).get(b).getSeqClusterId() + "   " + clusterList.get(a).get(b).getStrClusterId() + "   ");
		//				for(int c = 0; c < clusterList.get(a).get(b).getStrCluster().size(); c ++) {
		//					System.out.print(clusterList.get(a).get(b).getStrCluster().get(c)._1 + "   ");
		//				}
		//				System.out.println();
		//			}
		//		} 

		List<Tuple2<Integer, List<Integer>>> alteredSeqId = SequenceClusterIdComparison.getSeqIdClustering(sequenceIdentity, sc, startCluster, endCluster);
		List<Tuple2<Integer, List<Cluster>>> altSeqIdClusters = new ArrayList<Tuple2<Integer, List<Cluster>>>();
		for(int a = 0; a < alteredSeqId.size(); a ++) {
			altSeqIdClusters.add(new Tuple2<Integer, List<Cluster>> (alteredSeqId.get(a)._1, new ArrayList<Cluster>()));
			for(int b = 0; b < alteredSeqId.get(a)._2.size(); b++) {
				for(int c = 0; c < clusterList.size(); c++) {
					//					System.out.println(clusterList.get(c).get(0).getSeqClusterId() + ", " + alteredSeqId.get(a)._2.get(b));
					if(clusterList.get(c).get(0).getSeqClusterId() == alteredSeqId.get(a)._2.get(b).intValue()) {
						//						System.out.println(clusterList.get(c).get(0).getSeqClusterId() + ",,,,,, " + alteredSeqId.get(a)._2.get(b));
						altSeqIdClusters.get(a)._2.addAll(clusterList.get(c));
						break;
					}
				}
			}
		}

		for(int p = 0; p < altSeqIdClusters.size(); p ++) {
			for(int q = 0; q < altSeqIdClusters.get(p)._2.size(); q ++) {
				altSeqIdClusters.get(p)._2.get(q).findRepChain();
				altSeqIdClusters.get(p)._2.get(q).findStrRepChain();
			}
		}

		for(int x = 0; x < altSeqIdClusters.size(); x ++) {
			for(int y = 0; y < altSeqIdClusters.get(x)._2.size() - 1; y ++) {
				Cluster strCluster1 = altSeqIdClusters.get(x)._2.get(y);
				List<Tuple3<String, SimplePolymerChain, Double>> strRepChain1 = strCluster1.getStrRepChain();
				for(int z = y +1; z < altSeqIdClusters.get(x)._2.size(); z ++) {
					Cluster strCluster2 = altSeqIdClusters.get(x)._2.get(z);
					List<Tuple3<String, SimplePolymerChain, Double>> strRepChain2 = strCluster2.getStrRepChain();
					boolean canMerge = true;
					if(strCluster1.getSeqClusterId() == strCluster2.getSeqClusterId()) {
						canMerge = false;
					}
					else {
						for(int alpha = 0; alpha < strRepChain1.size(); alpha ++) {
							for(int beta = 0; beta < strRepChain2.size(); beta ++) {
								double rmsd = SmithWaterman.getFatCatTmScore(strRepChain1.get(alpha)._2().getCoordinates(), strRepChain2.get(beta)._2().getCoordinates())[1];
								if(strRepChain1.get(alpha)._3() + strRepChain2.get(beta)._3() + rmsd > maxRmsd) {
									//	System.out.println(strRepChain1.get(alpha)._3() + "  " + strRepChain2.get(beta)._3() + "  " + rmsd);
									canMerge = false;
								}
								else {
									//	System.out.println(strRepChain1.get(alpha)._3() + "  " + strRepChain2.get(beta)._3() + "  " + rmsd);
								}
							}
						}
					}
					if(canMerge == true) {
						Cluster mergedCluster = altSeqIdClusters.get(x)._2.get(y).merge(altSeqIdClusters.get(x)._2.get(z), 0, 0);
						altSeqIdClusters.get(x)._2.set(y, mergedCluster);
						altSeqIdClusters.get(x)._2.remove(z);
						z--;
					}
				}
			}
		}

		for(int j = 0; j < altSeqIdClusters.size(); j ++) {
			int tempSeqClusterId = altSeqIdClusters.get(j)._1;
			for(int k = 0; k < altSeqIdClusters.get(j)._2.size(); k ++) {
				altSeqIdClusters.get(j)._2.get(k).setSeqClusterId(tempSeqClusterId);
				altSeqIdClusters.get(j)._2.get(k).setStrClusterId(k + 1);
			}
		}

		quickSortAlt(altSeqIdClusters, 0, altSeqIdClusters.size() - 1);

		PrintWriter writer = new PrintWriter(outputFileName);
		writer.println("SequenceCluster,StructuralCluster,StructuralClusterSize,RepresentativeChain,OtherChains");
		for(int m = 0; m < altSeqIdClusters.size(); m ++) {
			for(int n = 0; n < altSeqIdClusters.get(m)._2.size(); n ++) {
				writer.print(altSeqIdClusters.get(m)._2.get(n).getSeqClusterId());
				writer.print(",");
				writer.print(altSeqIdClusters.get(m)._2.get(n).getStrClusterId());
				writer.print(",");
				int currentClusterSize = altSeqIdClusters.get(m)._2.get(n).size();
				writer.print(currentClusterSize);
				writer.print(",");
				String repChainName = altSeqIdClusters.get(m)._2.get(n).getRepChain()._1;
				writer.print(repChainName);
				List<Tuple2<String, SimplePolymerChain>> tempCluster = altSeqIdClusters.get(m)._2.get(n).getStrCluster(); 
				for(int x = 0; x < currentClusterSize; x ++) {
					Tuple2<String, SimplePolymerChain> tempChain = tempCluster.get(x);
					if(!tempChain._1.equals(repChainName)) {
						writer.print("," + tempCluster.get(x)._1);
					}
				}
				writer.println();
			}
		}
		writer.close();
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
				.sequenceFile(hadoopSequenceFileName, Text.class, SimplePolymerChain.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.mapToPair(new HadoopToSimplePolymerChainMapper())
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

	/**
	 * returns the total number of 100% sequence clusters in the PDB
	 * @return the number of sequence clusters
	 */
	private static int getNumSeqClusters() {
		int sequenceIdentity = 100;
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);

		List<List<String>> clusters = reader.getPdbChainIdClusters();
		return clusters.size();
	}

	//	/**
	//	 * returns the total number of sequence clusters with seqId identity in the PDB
	//	 * @param seqId the sequence identity 
	//	 * @return the number of sequence clusters
	//	 */
	//	private static int getNumSeqClusters(int seqId) {
	//		int sequenceIdentity = seqId;
	//		BlastClustReader reader = new BlastClustReader(sequenceIdentity);
	//
	//		List<List<String>> clusters = reader.getPdbChainIdClusters();
	//		return clusters.size();
	//	}

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

	/**
	 * Uses quickSort algorithm to sort a list of clusters based on their sequence cluster Ids
	 * @param strClusters the list of clusters to be sorted
	 * @param first the index of the first cluster in the list of clusters to be sorted
	 * @param last the index of the last cluster in the list of clusters to be sorted
	 */
	public void quickSort(List<Cluster> strClusters, int first, int last) {
		int g = first, h = last;
		int midIndex = (first + last) / 2;
		int dividingValue = strClusters.get(midIndex).getSeqClusterId();
		int dividingValue2 = strClusters.get(midIndex).getStrClusterId();
		do {
			while((strClusters.get(g).getSeqClusterId() < dividingValue) || (strClusters.get(g).getSeqClusterId() == dividingValue && strClusters.get(g).getStrClusterId() < dividingValue2)) g++;
			while((strClusters.get(h).getSeqClusterId() > dividingValue) || (strClusters.get(h).getSeqClusterId() == dividingValue && strClusters.get(h).getStrClusterId() > dividingValue2)) h--;
			if(g <= h) {
				Cluster temp = strClusters.get(g);
				strClusters.set(g, strClusters.get(h));
				strClusters.set(h, temp);
				g++;
				h--;
			}
		}
		while (g < h);
		if (h > first) quickSort (strClusters, first, h);
		if (g < last) quickSort (strClusters, g, last);
	}

	/**
	 * Uses quickSort algorithm to sort a List<Tuple2<Integer, List<Cluster>>> based on their sequence cluster Ids
	 * @param altSeqIdClusters the List<Tuple2<Integer, List<Cluster>>> to be sorted
	 * @param first the index of the first Tuple2 to be sorted
	 * @param last the index of the last Tuple2 to be sorted
	 */
	public void quickSortAlt(List<Tuple2<Integer, List<Cluster>>> altSeqIdClusters, int first, int last) {
		int g = first, h = last;
		int midIndex = (first + last) / 2;
		int dividingValue = altSeqIdClusters.get(midIndex)._1;
		do {
			while(altSeqIdClusters.get(g)._1 < dividingValue) g++;
			while(altSeqIdClusters.get(h)._1 > dividingValue) h--;
			if(g <= h) {
				Tuple2<Integer,List<Cluster>> temp = altSeqIdClusters.get(g);
				altSeqIdClusters.set(g, altSeqIdClusters.get(h));
				altSeqIdClusters.set(h, temp);
				g++;
				h--;
			}
		}
		while (g < h);
		if (h > first) quickSortAlt (altSeqIdClusters, first, h);
		if (g < last) quickSortAlt (altSeqIdClusters, g, last);
	}

	//	/**
	//	 * Uses quickSort algorithm to sort a list of clusters based on their structural cluster Ids
	//	 * @param strClusters the list of clusters to be sorted
	//	 * @param first the index of the first cluster in the list of clusters to be sorted
	//	 * @param last the index of the last cluster in the list of clusters to be sorted
	//	 */
	//	public void quickSortByStrId(List<Cluster> strClusters, int first, int last) {
	//		int g = first, h = last;
	//		int midIndex = (first + last) / 2;
	//		int dividingValue = strClusters.get(midIndex).getStrClusterId();
	//		do {
	//			while(strClusters.get(g).getStrClusterId() < dividingValue) g++;
	//			while(strClusters.get(h).getStrClusterId() > dividingValue) h--;
	//			if(g <= h) {
	//				Cluster temp = strClusters.get(g);
	//				strClusters.set(g, strClusters.get(h));
	//				strClusters.set(h, temp);
	//				g++;
	//				h--;
	//			}
	//		}
	//		while (g < h);
	//		if (h > first) quickSort (strClusters, first, h);
	//		if (g < last) quickSort (strClusters, g, last);
	//	}
}
