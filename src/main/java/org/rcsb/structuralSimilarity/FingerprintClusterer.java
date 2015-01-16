package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.clustering.KMeans;
import org.apache.spark.mllib.clustering.KMeansModel;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.rcsb.fingerprints.DCT1DFingerprint;

import scala.Tuple2;

/**
 * 
 * @author  Peter Rose
 */
public class FingerprintClusterer 
{    
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 4;
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = "/Users/peter/Data/PDB_CHAINS/protein_chains_40_20150111_232623.seq";
		String results = "/Users/peter/Data/ProteinSimilarity/FingerPrintClusterer20150104.csv";

		long t1 = System.nanoTime();

		SparkConf conf = new SparkConf().setMaster("local[*]").setAppName("FingerprintClusterer");
		conf.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
		
		// calculate fingerprint vectors
		JavaPairRDD<String, Vector> chains = sc.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS)  
//				.sample(false, 0.1, 123456)
		        .mapToPair(new SeqToChainMapper()) // convert input chain id, CA coordinate pairs
		        .mapToPair(new ChainToFeatureVectorMapper(new DCT1DFingerprint()))  // calculate fingerprints
		        .cache();

		List<Tuple2<String,Vector>> bc = chains.collect();
		final Broadcast<List<Tuple2<String,Vector>>> vec = sc.broadcast(bc);
		
		JavaRDD<Vector> vectors = chains.values().cache();
		JavaRDD<Vector> vsqrt = vectors.map(new Function<Vector, Vector>() {
			private static final long serialVersionUID = 1L;

			@Override
			public Vector call(Vector v) throws Exception {
				double[] array = v.toArray();
				for (int i = 0; i < array.length; i++) {
					array[i] = Math.sqrt(array[i]);
				}
				// TODO Auto-generated method stub
				return Vectors.dense(array);
			}
		});
		chains.unpersist();
		
		int numVectors = bc.size();
		System.out.println("Vectors: " + numVectors);
		
		long t2 = System.nanoTime();

		// KMeans clustering by fingerprint vectors		
		int numIterations = 20;
		// assume less than 5% of total are similar
		int numClusters = (int)Math.round(numVectors*0.1);
		System.out.println("numClusters: " + numClusters);

		KMeans km = new KMeans();
		km.setInitializationMode(KMeans.K_MEANS_PARALLEL());
//		KMeansModel clusters = KMeans.train(vectors.rdd(), numClusters, numIterations);
//		double wssse = clusters.computeCost(vectors.rdd());
		KMeansModel clusters = KMeans.train(vsqrt.rdd(), numClusters, numIterations);
		double wssse = clusters.computeCost(vsqrt.rdd());
		System.out.println("clusters: " + numClusters + " error: " + wssse);

//		JavaRDD<Integer> membership = clusters.predict(vectors);
		JavaRDD<Integer> membership = clusters.predict(vsqrt);
		
		long t3 = System.nanoTime();
		
		// create pairs within each cluster
		List<Tuple2<Integer, Integer>> pairList = createPairList(membership);
		long numPairs = pairList.size();
		System.out.println("number of pairs: " + numPairs);
		
		List<Tuple2<String, Float>> bestScores = sc.parallelizePairs(pairList, NUM_THREADS*NUM_TASKS_PER_THREAD)
				.mapToPair(new PairSimilarityCalculator(vec))
			    .filter(s -> s._2 > 0.9f)
			    .collect();
		
		long t4 = System.nanoTime();  
		
		// save results
		writeToCsv(results, bestScores);
		sc.stop();
		sc.close();

		System.out.println("numClusters: " + numClusters);
		System.out.println(" error: " + wssse);
		System.out.println("pairs: " + numPairs);
		System.out.println("filtered pairs: " + bestScores.size());
		System.out.println("calculate fingerprints : " + (t2 - t1)/1E9 + " s");
		System.out.println("KMeans clustering      : " + (t3-t2)/1E6 + " s");
		System.out.println("calculate pairs        : " + (t4-t3)/1E6 + " s");
		System.out.println("Time: " + (t4 - t1)/1E9 + " sec.");
	}

	private static void writeToCsv(String results, List<Tuple2<String, Float>> bestScores) throws FileNotFoundException {
		PrintWriter writer = new PrintWriter(results);
		for (Tuple2<String, Float> t : bestScores) {
			writer.print(t._1);
			writer.print(",");
			writer.print(t._2);
			writer.println();
		}
		writer.flush();
		writer.close();
	}

	private static List<Tuple2<Integer, Integer>> createPairList(JavaRDD<Integer> membership) {
		JavaPairRDD<Integer, Long> mm = membership.zipWithIndex();		
		JavaPairRDD<Integer, Iterable<Long>> group = mm.groupByKey();
		
		List<Tuple2<Integer, Integer>> pairList = new ArrayList<>();
		for (Tuple2<Integer,Iterable<Long>> cluster: group.collect()) {	
			List<Tuple2<Integer, Integer>> list = toPairs(cluster._2);
			if (! list.isEmpty()) {
				pairList.addAll(list);
			}
		}
		return pairList;
	}
	
	private static List<Tuple2<Integer, Integer>> toPairs(Iterable<Long> cluster) {
		List<Integer> members = new ArrayList<>();
		for (Long member: cluster) {
			members.add(member.intValue());
		}
		int n = members.size();
		List<Tuple2<Integer,Integer>> list = new ArrayList<>((n*(n-1)/2));
		for (int i = 0; i < n-1; i++) {
			for (int j = i+1; j < n; j++) {
				list.add(new Tuple2<Integer,Integer>(members.get(i),members.get(j)));
			}
		}
		return list;
	}
	
//	class MyRegistrator extends KryoRegistrator {
//		  override def registerClasses(kryo: Kryo) {
//		    kryo.register(classOf[MyClass1])
//		    kryo.register(classOf[MyClass2])
//		  }
//		}

}

