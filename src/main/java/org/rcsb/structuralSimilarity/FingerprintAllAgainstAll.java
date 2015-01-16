package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.rdd.HadoopRDD;
import org.rcsb.fingerprints.DCT1DFingerprint;
import org.rcsb.fingerprints.EndToEndDistanceFingerprint;

import scala.Tuple2;

/**
 * 
 * @author  Peter Rose
 */
public class FingerprintAllAgainstAll {  
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 4;
	private static int BATCH_SIZE = 5000000;

	private int pairsProcessed = 0;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = "/Users/peter/Data/PDB_CHAINS/protein_chains_40_20150114_141156.seq";
		String results = "/Users/peter/Data/ProteinSimilarity/FingerPrintAllAgainstAll20150111.csv";

		FingerprintAllAgainstAll aaa = new FingerprintAllAgainstAll();
		aaa.run(path, results);
	}

	private void run(String path, String results) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName("1" + this.getClass().getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer").setSparkHome("/tmp/");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		long t1 = System.nanoTime();

		// read file with protein chains and calculate fingerprint vectors
		List<Tuple2<String,Vector>> bc  = sc.
				sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS)  // read protein chains
//				.sample(false, 0.4, 123456)
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate> pairs
				.filter(new GapFilter(3, 5))
				.filter(new LengthFilter(75,1000))
				.mapToPair(new ChainToFeatureVectorMapper(new EndToEndDistanceFingerprint())) // calculate fingerprints
	//			.mapToPair(new ChainToFeatureVectorMapper(new DCT1DFingerprint())) // calculate fingerprints
				.collect();

		// broadcast feature vectors to all nodes
		final Broadcast<List<Tuple2<String,Vector>>> vec = sc.broadcast(bc);
		int numVectors = bc.size();
		System.out.println("Vectors: " + numVectors);

		long t2 = System.nanoTime();

		// process pairwise comparisons in batches
		System.out.println("Total pairs: " + (numVectors*(numVectors-1)/2));
		
		PrintWriter writer = new PrintWriter(results);	

		int numBestScores = 0;
		long numPairs = 0;
		this.pairsProcessed = 0;
		
		while (this.pairsProcessed < numVectors) {
			System.out.println("Pairs processed: " + numPairs);

			// get a batch of pairs
			List<Tuple2<Integer, Integer>> pairList = nextBatch(numVectors);

			// calculate pairwise scores and filter with score > 0.9
			List<Tuple2<String, Float>> list = sc.parallelizePairs(pairList, NUM_THREADS*NUM_TASKS_PER_THREAD)
					.mapToPair(new PairSimilarityCalculator(vec))
					.filter(s -> s._2 > 0.9f)
					.collect();	

			numPairs += pairList.size();
			numBestScores += list.size();
			// write results to .csv file
			writeToCsv(writer, list);
		}
		long t3 = System.nanoTime();

		writer.close();
		sc.stop();
		sc.close();

		System.out.println("protein chains: " + numVectors);
		System.out.println("pairs         : " + numPairs);
		System.out.println("filtered pairs: " + numBestScores);
		System.out.println();
		System.out.println("calculate fingerprints : " + (t2-t1)/1E9 + " s");
		System.out.println("calculate pairs        : " + (t3-t2)/1E9 + " s");
		System.out.println("total time             : " + (t3 - t1)/1E9 + " s");
		System.out.println("time per pair          : " + ((t3 - t1)/numPairs)  + " ns");
	}

	private static void writeToCsv(PrintWriter writer, List<Tuple2<String, Float>> list) {
		for (Tuple2<String, Float> t : list) {
			writer.print(t._1);
			writer.print(",");
			writer.print(t._2);
			writer.println();
		}
		writer.flush();
	}

	private List<Tuple2<Integer, Integer>> nextBatch(int n) {
		List<Tuple2<Integer,Integer>> list = new ArrayList<>(BATCH_SIZE);

		for (int i = this.pairsProcessed; i < n-1; i++) {
			for (int j = i+1; j < n; j++) {
				list.add(new Tuple2<Integer,Integer>(i,j));
			}
			if (list.size() > BATCH_SIZE*0.9) {
				this.pairsProcessed = i + 1;
				return list;
			}
		}
		this.pairsProcessed = n;
		return list;
	}
}

