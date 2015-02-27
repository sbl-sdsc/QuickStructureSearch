package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;
import org.rcsb.fingerprints.CombinationFingerprint;
import org.rcsb.fingerprints.DCT1DFingerprint;
import org.rcsb.fingerprints.DCT1DLinearFingerprint;
import org.rcsb.fingerprints.DCT1DOptFingerprint;
import org.rcsb.fingerprints.EndToEndDistanceFingerprint;
import org.rcsb.fingerprints.GenericFingerprint;
import org.rcsb.fingerprints.LinearFingerprint;
import org.rcsb.fingerprints.NullHypothesisFingerprint;
import org.rcsb.fingerprints.PointToPointDistanceFingerprint;
import org.rcsb.fingerprints.TetrahedronFingerprint;

import scala.Tuple2;

/**
 * This class reads .csv files of protein chain pairs and TM similarity metrics and compares
 * the results with similarities calculated by fingerprinting methods
 * 
 * @author  Peter Rose
 */
public class DCT1DOpt_Optimizer { 
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		if (args.length < 2) {
			System.out.println("Usage: FingerPrintTester.jar inputDirectory [sequenceFile] outputFile");
			System.out.println("  inputDirectory: directory with .csv files that contain TM scores");
			System.out.println("  sequenceFile: Hadoop sequence file with protein chain information");
			System.out.println("  outputFile: results from calculation in .csv format. This file should have a .csv extension");
		}
	
		String sequenceFileName = args[1]; 
	    String outputFileName = args[2];
	    int randomSeed = Integer.parseInt(args[3]);
	
		DCT1DOpt_Optimizer aaa = new DCT1DOpt_Optimizer();
		aaa.run(args[0], sequenceFileName, outputFileName, randomSeed);
	}

	private void run(String inputDirectories, String path, String outputFileName, int randomSeed) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		PrintWriter writer = new PrintWriter(outputFileName);
		
		long t1 = System.nanoTime();
		
		// split input lines of .csv files into chainId1,chainId2, data pairs
        JavaPairRDD<String, String> trainingData = sc
        		.textFile(inputDirectories, NUM_THREADS*NUM_TASKS_PER_THREAD)
        		.mapToPair(new SplitAtComma(2)).cache();
         
        // create list of chain id pairs <chainId1,chainId2>
        JavaPairRDD<String, String> pairs = trainingData
           		.keys()
           		.mapToPair(new SplitAtComma(1))
           		.cache();
        
        // broadcast list of unique chain ids
        List<String> chainIds = pairs
        		.keys().distinct() // all distinct keys (chainId1)
        		.union(pairs.values().distinct()) // union with all distinct values (chainId2)
        		.distinct() // finally make sure all chain ids are distinct
        		.collect();
        final Broadcast<Set<String>> chainIdsBc = sc.broadcast(new HashSet<String>(chainIds));



        // calculate <chainId, feature vector> pairs
        JavaPairRDD<String, Point3d[]> cache = sc
        		.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS*NUM_TASKS_PER_THREAD)  // read protein chains
        		//	.sample(false, 0.1, 123456) // use only a random fraction, i.e., 40%
        		.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate[]> pairs				
        		.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 3 and <= 5 gaps
        		.filter(new LengthFilter(50,500)) // keep protein chains with at least 50 residues
        		.filter(new ChainIdFilter<Point3d[]>(chainIdsBc)) // calculate feature vectors for chains in the training set only
        		//		     	.mapToPair(new ChainSmootherMapper(new RogenChainSmoother(2))) // add new chain smoother here ...
        		.mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(1))) // add new chain smoother here ...
        		.cache();

        int nTrials = 300; // 19170 sec.
        Random r = new Random(randomSeed);
        
        for (int trial = 0; trial < nTrials; trial++) {

//        	DCT1DOptFingerprint fp = new DCT1DOptFingerprint(r.nextInt());
        	DCT1DLinearFingerprint fp = new DCT1DLinearFingerprint(r.nextInt());
        	float[] parameters = fp.getParameters();
        	JavaPairRDD<String, Vector> features = cache
//        			.mapToPair(new ChainToFeatureVectorMapper(fp))
        		    .mapToPair(new ChainToLinearFeatureMapper(fp))
        			.cache();// calculate features

        	// broadcast feature vectors
        	List<Tuple2<String,Vector>>  bc =  features.collect(); // return results to master node     
        	final Broadcast<List<Tuple2<String,Vector>>> featureVectorsBc = sc.broadcast(bc);

        	// broadcast list of chain ids that have feature vectors
        	List<String> availableChainIds = features
        			.keys()
        			.collect();
        	final Broadcast<List<String>> availableChainIdsBc = sc.broadcast(availableChainIds);

        	// calculate Jaccard Index and join with TM metrics
        	List<Tuple2<String, Tuple2<Float, String>>> results = pairs
        			.filter(new ChainIdPairFilter(availableChainIdsBc)) // only keep pairs that have feature vectors available
        			.mapToPair(new ChainIdToIndexMapper(availableChainIdsBc)) // map chain ids to indices into feature vector
 //       			.mapToPair(new FeatureVectorToJaccardMapper(featureVectorsBc)) // maps pairs of feature vectors to Jaccard index
        			.mapToPair(new LinearFeatureVectorToLevenshteinMapper(featureVectorsBc))
        			//			.mapToPair(new FeatureVectorToContainmentScoreMapper(featureVectorsBc)) // maps pairs of feature vectors to Jaccard index
        			.join(trainingData) // join with TM metrics from the input file
        			.sortByKey()
        			.collect();

        	float[] scores = getStatistics(results, 0.5f , 0.5f);

        	writeToCsv2(writer, parameters, scores);
        	System.out.println("Trial: " + trial);
        	System.out.println(parameters + " --- " + scores);
        	System.out.println("--------------------------------------");
        	features.unpersist();
        }
        sc.stop();
        sc.close();

        writer.close();


        long t3 = System.nanoTime();
        System.out.println();

        System.out.println("total time         : " + (t3-t1)/1E9 + " s");
	}

	/**
	 * Writes parameters and statistics to a csv file
	 * @param writer
	 * @param list
	 */
	private static void writeToCsv2(PrintWriter writer, float[] parameters, float[] statistics) {
		for (float f: parameters) {
			writer.print(f); // chainId pair
			writer.print(",");
		}
		for (int i = 0; i < statistics.length-1; i++) {
			writer.print(statistics[i]); // chainId pair
			writer.print(",");
		}
		writer.print(statistics[statistics.length-1]);
		writer.println();
		writer.flush();
	}

	/**
	 * 
	 * @param writer
	 * @param list
	 */
	private static void printStatistics(List<Tuple2<String, Tuple2<Float, String>>> joinedResults) {	
		System.out.printf("%7s %7s %7s %7s %7s %7s %7s %7s", "F", "TP", "FN", "TN", "FP", "SENS", "SPEC", "F1");
		System.out.println();
		float tmThreshold = 0.5f;
		for (float f = 0.3f; f < 0.8f; f+= 0.05f) {
			float[] scores = getStatistics(joinedResults, f, tmThreshold);
			System.out.printf("%8.2f", f);
			System.out.printf("%8d", (int)scores[0]);
			System.out.printf("%8d", (int)scores[1]);
			System.out.printf("%8d", (int)scores[2]);
			System.out.printf("%8d", (int)scores[3]);
			System.out.printf("%8.2f", scores[4]);
			System.out.printf("%8.2f", scores[5]);
			System.out.printf("%8.2f", scores[6]);
			System.out.println();
			//			System.out.println(f + ": " + Arrays.toString(scores));
		}
	}

	private static float[] getStatistics(List<Tuple2<String, Tuple2<Float, String>>> joinedResults, float threshold, float tmThreshold) {
		float[] scores = new float[7];

		int tp = 0;
		int tn = 0;
		int fp = 0;
		int fn = 0;

		for (Tuple2<String, Tuple2<Float, String>> t: joinedResults) {
			float tmScore = Float.parseFloat(t._2._2.split(",")[0]);
			float fingerPrintScore = t._2._1;
			if (tmScore >= tmThreshold) {
				if (fingerPrintScore >= threshold) {
					tp++;
				} else {
					fn++;
				}
			} else {
				if (fingerPrintScore >= threshold)	 {
					fp++;
				} else {
					tn++;
				}
			}
		}
		scores[0] = tp;
		scores[1] = fn;
		scores[2] = tn;
		scores[3] = fp;
		scores[4] = tp/(float)(tp+fn);
		scores[5] = tn/(float)(fp+tn);
		scores[6] = 2*tp/(float)(2*tp+fp+fn);

		return scores;
	}
}

