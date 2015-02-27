package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;
import org.rcsb.fingerprints.DCT1DFingerprint;
import org.rcsb.fingerprints.EndToEndDistanceFingerprint;
import org.rcsb.fingerprints.TetrahedronFingerprint;

import scala.Tuple2;

/**
 * 
 * @author  Peter Rose
 */
public class FingerprintAllAgainstAll { 
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 5000000; // number of pairs to be processed per batch

	private int pairsProcessed;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = "src/test/resources/protein_chains_40_20150114_141156.seq";
		String outputFileName = args[0];

		if (args.length == 2) {
			sequenceFileName = args[0]; 
			outputFileName = args[1];
		}
		FingerprintAllAgainstAll aaa = new FingerprintAllAgainstAll();
		aaa.run(sequenceFileName, outputFileName);
	}

	private void run(String path, String results) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		long t1 = System.nanoTime();

		// Step 1. calculate <pdbId.chainId, feature vector> pairs
		List<Tuple2<String,Vector>> bc  = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS)  // read protein chains
			//	.sample(false, 0.1, 123456) // use only a random fraction, i.e., 40%
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate> pairs
				.filter(new GapFilter(0, 5)) // keep protein chains with gap size <= 3 and <= 5 gaps
				.filter(new LengthFilter(50,1000)) // keep protein chains with at least 50 residues
		   //	.mapToPair(new ChainSmootherMapper(new RogenChainSmoother(2))); // add new chain smoother here ...
				.mapToPair(new ChainToFeatureVectorMapper(new TetrahedronFingerprint())) // calculate features
			//	.mapToPair(new ChainToFeatureVectorMapper(new EndToEndDistanceFingerprint())) // calculate features
	       	//  .mapToPair(new ChainToFeatureVectorMapper(new DCT1DFingerprint())) // calculate features
				.collect(); // return results to master node

		// Step 2.  broadcast feature vectors to all nodes
		final Broadcast<List<Tuple2<String,Vector>>> featureVectors = sc.broadcast(bc);
		int numVectors = bc.size();
		System.out.println("Vectors: " + numVectors);

		long t2 = System.nanoTime();

		// Step 3: process pairwise comparisons in batches
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
			List<Tuple2<String, Float>> list = sc
					.parallelizePairs(pairList, NUM_THREADS*NUM_TASKS_PER_THREAD) // distribute data
					.mapToPair(new FeatureVectorToJaccardMapper(featureVectors)) // maps pairs of feature vectors to Jaccard index
					.filter(s -> s._2 > 0.9f) // keep only a pair with a Jaccard index > 0.9
					.collect();	// copy result to master node

			// write results to .csv file
			writeToCsv(writer, list);
			
			numPairs += pairList.size();
			numBestScores += list.size();
		}
		long t3 = System.nanoTime();

		writer.close();
		sc.stop();
		sc.close();

		System.out.println("protein chains     : " + numVectors);
		System.out.println("total pairs        : " + numPairs);
		System.out.println("filtered pairs     : " + numBestScores);
		System.out.println();
		System.out.println("calculate features : " + (t2-t1)/1E9 + " s");
		System.out.println("compare pairs      : " + (t3-t2)/1E9 + " s");
		System.out.println("total time         : " + (t3-t1)/1E9 + " s");
		System.out.println("time per pair      : " + ((t3-t1)/numPairs)  + " ns");
	}

	/**
	 * Writes pairs of chain ids and calculated similarity score to a csv file
	 * @param writer
	 * @param list
	 */
	private static void writeToCsv(PrintWriter writer, List<Tuple2<String, Float>> list) {
		for (Tuple2<String, Float> t : list) {
			writer.print(t._1);
			writer.print(",");
			writer.print(t._2);
			writer.println();
		}
		writer.flush();
	}

	/**
	 * Returns pairs of indices for the pairwise comparison. This is done 
	 * in batches to reduce the memory footprint.
	 * @param n number of feature vectors
	 * @return
	 */
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

