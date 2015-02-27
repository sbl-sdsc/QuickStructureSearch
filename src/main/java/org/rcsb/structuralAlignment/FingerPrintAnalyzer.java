package org.rcsb.structuralAlignment;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaDoubleRDD;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.classification.NaiveBayes;
import org.apache.spark.mllib.classification.NaiveBayesModel;
import org.apache.spark.mllib.evaluation.BinaryClassificationMetrics;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.regression.LabeledPoint;
import org.apache.spark.mllib.stat.MultivariateStatisticalSummary;
import org.apache.spark.mllib.stat.Statistics;
import org.rcsb.fingerprints.DCT1DFingerprint;
import org.rcsb.fingerprints.DCT1DOptFingerprint;
import org.rcsb.fingerprints.EndToEndDistanceFingerprint;
import org.rcsb.fingerprints.GenericFingerprint;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;
import org.rcsb.structuralSimilarity.SeqToChainMapper;

import scala.Tuple2;

/**
 * This class creates structural alignments between random protein chain pairs 
 * using jFatCAT and scores the alignments with the TM score
 * 
 * @author  Peter Rose
 */
public class FingerPrintAnalyzer implements Serializable { 
	private static final long serialVersionUID = 2779213801755875110L;
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 100000;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0]; 
		String outputFileName = args[1];
		int nPairs = Integer.parseInt(args[2]);
		int seed = Integer.parseInt(args[3]);
		
		long t1 = System.nanoTime();
		FingerPrintAnalyzer creator = new FingerPrintAnalyzer();
		creator.run(sequenceFileName, outputFileName, nPairs, seed);
		System.out.println("Time: " + ((System.nanoTime()-t1)/1E9) + " s");
	}

	private void run(String path, String outputFileName, int nPairs, int seed) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g");
//				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		// Step 1. calculate <pdbId.chainId, feature vector> pairs
        List<Tuple2<String, Point3d[]>> chains = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS)  // read protein chains
				.sample(false, 0.1, 123456) // use only a random fraction, i.e., 40%
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate> pairs
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 3 and <= 5 gaps
				.filter(new LengthFilter(50,500)) // keep protein chains with at least 50 residues
//				.mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(1))) // add new chain smoother here ...
				.collect(); // return results to master node

		// Step 2.  broadcast feature vectors to all nodes
		final Broadcast<List<Tuple2<String,Point3d[]>>> chainsBc = sc.broadcast(chains);
		int nChains = chains.size();

		Random r = new Random(seed);

		PrintWriter writer = new PrintWriter(outputFileName);
//		GenericFingerprint fingerprint = new DCT1DFingerprint(8, 40);
		int fragmentLength = 9;
		GenericFingerprint fingerprint = new DCT1DOptFingerprint();
//		GenericFingerprint fingerprint = new EndToEndDistanceFingerprint(8, 1);
		
		List<Tuple2<Integer, Integer>> pairs = randomPairs(nChains, nPairs, r.nextLong());

		JavaRDD<LabeledPoint> data = sc
				.parallelizePairs(pairs, NUM_THREADS * NUM_TASKS_PER_THREAD)
				.map(new ChainPairToRmsdMapper1(chainsBc, fingerprint, fragmentLength))
				.cache(); // copy result to master node

	
		JavaRDD<LabeledPoint> positives = data.filter(p -> p.label() > 0).cache();
		JavaRDD<LabeledPoint> negatives = data.filter(p -> p.label() < 1).cache();
		
		JavaRDD<Vector> sPos = positives.map(p -> p.features());
		MultivariateStatisticalSummary summary = Statistics.colStats(sPos.rdd());
		double pMean = summary.mean().toArray()[0]; // a dense vector containing the mean value for each column
		double pVariance = summary.variance().toArray()[0]; // column-wise variance	
		
		JavaRDD<Vector> sNeg = negatives.map(p -> p.features());
		summary = Statistics.colStats(sNeg.rdd());
		double nMean = summary.mean().toArray()[0]; // a dense vector containing the mean value for each column
		double nVariance = summary.variance().toArray()[0]; // column-wise variance
		
		
		
		long pos = positives.count();
		long neg = negatives.count();	
		
	    JavaRDD<Double> pValues = positives.map(p -> p.features().toArray()[0]);
		Double maxValue = pValues.max(Comparator.naturalOrder());
		
		long mis = negatives.filter(p -> p.features().toArray()[0] < maxValue).count();
				
		writeToCsv(writer, data.collect());
		writer.close();
		
		sc.stop();
		sc.close();

		System.out.println("protein chains     : " + nChains);
		System.out.println("ramdom pairs        : " + nPairs);
		System.out.println("maxvalue  : " + maxValue);
		System.out.println("positives : " + pos);
	    System.out.println("negatives : " + neg);
	    System.out.println("pMean     : " + pMean + " +- " + pVariance);
	    System.out.println("nMean     : " + nMean + " +- " + nVariance);
	    System.out.println("mismatches: " + mis);
	}

	/**
	 * Writes pairs to a csv file
	 * @param list
	 * @param list
	 */
	private static void writeToCsv(PrintWriter writer, List<LabeledPoint> list) {
		for (LabeledPoint p : list) {
			writer.print(p.label());
			writer.print(",");
			writer.println(p.features().toArray()[0]);
		}
		// writer.flush();
	}

	/**
	 * Returns random pairs of indices for the pairwise comparison.
	 * @param n number of feature vectors
	 * @return
	 */
	private List<Tuple2<Integer, Integer>> randomPairs(int n, int nPairs, long seed) {
		Random r = new Random(seed);
		List<Tuple2<Integer,Integer>> list = new ArrayList<>(nPairs);

		for (int i = 0; i < nPairs; i++) {
			int j = r.nextInt(n);
			int k = r.nextInt(n);
			if (j == k) {
				continue;
			}

			Tuple2<Integer,Integer> tuple = new Tuple2<>(j,k);
			if (! list.contains(tuple)) {
			    list.add(tuple);
			}
		}
		return list;
	}
}

