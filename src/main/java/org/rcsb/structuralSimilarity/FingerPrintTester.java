package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;
import org.rcsb.fingerprints.DCT1DFingerprint;
import org.rcsb.fingerprints.EndToEndDistanceFingerprint;
import org.rcsb.fingerprints.TetrahedronFingerprint;

import scala.Tuple2;

/**
 * This class reads .csv files of protein chain pairs and TM similarity metrics and compares
 * the results with similarities calculated by fingerprinting methods
 * 
 * @author  Peter Rose
 */
public class FingerPrintTester { 
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		if (args.length < 2) {
			System.out.println("Usage: CompareWithTm.jar inputDirectory [sequenceFile] outputFile");
			System.out.println("  inputDirectory: directory with .csv files that contain TM scores");
			System.out.println("  sequenceFile: Hadoop sequence file with protein chain information");
			System.out.println("  outputFile: results from calculation in .csv format. This file should have a .csv extension");
		}
		String sequenceFileName = "src/test/resources/protein_chains_40_20150114_141156.seq";
		String outputFileName = args[1];

		if (args.length == 3) {
			sequenceFileName = args[1]; 
			outputFileName = args[2];
		}
		FingerPrintTester aaa = new FingerPrintTester();
		aaa.run(args[0], sequenceFileName, outputFileName);
	}

	private void run(String inputDirectories, String path, String outputFileName) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName("1" + this.getClass().getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
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
        final Broadcast<List<String>> chainIdsBc = sc.broadcast(chainIds);
        
		// calculate <chainId, feature vector> pairs
        JavaPairRDD<String, Vector> features = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS*NUM_TASKS_PER_THREAD)  // read protein chains
			//	.sample(false, 0.1, 123456) // use only a random fraction, i.e., 40%
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate[]> pairs
				.filter(new ChainIdFilter<Point3d[]>(chainIdsBc)) // calculate feature vectors for chains in the training set only
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 3 and <= 5 gaps
				.filter(new LengthFilter(50,1000)) // keep protein chains with at least 50 residues
//		     	.mapToPair(new ChainSmootherMapper(new RogenChainSmoother(2))); // add new chain smoother here ...
//				.mapToPair(new ChainToFeatureVectorMapper(new TetrahedronFingerprint())) // calculate features
//				.mapToPair(new ChainToFeatureVectorMapper(new EndToEndDistanceFingerprint())) // calculate features
	       	    .mapToPair(new ChainToFeatureVectorMapper(new DCT1DFingerprint())) // calculate features
				.cache();
      
        // broadcast feature vectors
        List<Tuple2<String,Vector>>  bc =  features.collect(); // return results to master node     
		final Broadcast<List<Tuple2<String,Vector>>> featureVectorsBc = sc.broadcast(bc);
        int numVectors = bc.size();
		
		// broadcast list of chain ids that have feature vectors
		List<String> availableChainIds = features
				.keys()
				.collect();
		final Broadcast<List<String>> availableChainIdsBc = sc.broadcast(availableChainIds);
		
		long t2 = System.nanoTime();
		
		// calculate Jaccard Index and join with TM metrics
	    List<Tuple2<String, Tuple2<Float, String>>> results = pairs
				.filter(new ChainIdPairFilter(availableChainIdsBc)) // only keep pairs that have feature vectors available
				.mapToPair(new ChainIdToIndexMapper(availableChainIdsBc)) // map chain ids to indices into feature vector
				.mapToPair(new FeatureVectorToJaccardMapper(featureVectorsBc)) // maps pairs of feature vectors to Jaccard index
	//			.filter(s -> s._2 > 0.9f) // keep only a pair with a Jaccard index > 0.9
				.join(trainingData) // join with TM metrics from the input file
				.sortByKey()
				.collect();
		
		
		int numPairs = results.size();
		
		// write results to .csv file
		PrintWriter writer = new PrintWriter(outputFileName);
		writeToCsv2(writer, results);
		writer.close();

		long t3 = System.nanoTime();

		sc.stop();
		sc.close();
     	System.out.println("protein chains     : " + numVectors);
		System.out.println("total pairs        : " + numPairs);
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
	private static void writeToCsv2(PrintWriter writer, List<Tuple2<String, Tuple2<Float, String>>> joinedResults) {
		for (Tuple2<String, Tuple2<Float, String>> t: joinedResults) {
			writer.print(t._1); // chainId pair
			writer.print(",");
			writer.print(t._2._1); // fingerprint score
			writer.print(",");
			writer.print(t._2._2); // tm score and other related metrics
			writer.println();
		}
		writer.flush();
	}
}

