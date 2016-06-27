package org.rcsb.project3;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.hadoop.io.HadoopToSimplePolymerChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralSimilarity.ChainIdFilter;
import org.rcsb.structuralSimilarity.ChainIdPairFilter;
import org.rcsb.structuralSimilarity.ChainIdToIndexMapper;
import org.rcsb.structuralSimilarity.ChainSmootherMapper;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;
import org.rcsb.structuralSimilarity.SplitAtComma;

import scala.Tuple2;

/**
 * This class reads .csv files of protein chain pairs and TM similarity metrics and compares
 * the results one .csv file with similarities calculated by fingerprint methods
 * 
 * It is used for fingerprint and alignment testing
 * 
 * @author  Peter Rose
 */
public class FingerPrintTesterWithHadoop { 
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
		String sequenceFileName = "src/test/resources/protein_chains_All_20150629_002251.seq";
		String outputFileName = args[1];

		if (args.length == 3) {
			sequenceFileName = args[1]; 
			outputFileName = args[2];
		}
		FingerPrintTesterWithHadoop aaa = new FingerPrintTesterWithHadoop();
		aaa.run(args[0], sequenceFileName, outputFileName);
	}

	private void run(String inputDirectories, String path, String outputFileName) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
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
        final Broadcast<Set<String>> chainIdsBc = sc.broadcast(new HashSet<String>(chainIds));
        
		// calculate <chainId, feature vector> pairs
        JavaPairRDD<String, SequenceFeatureInterface<?>> features = sc
				.sequenceFile(path, Text.class, SimplePolymerChain.class, NUM_THREADS*NUM_TASKS_PER_THREAD)  // read protein chains
				.mapToPair(new HadoopToSimplePolymerChainMapper()) // convert input to <pdbId.chainId, protein sequence> pairs
				.filter(t -> t._2.isProtein())
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))				
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 3 and <= 5 gaps
				.filter(new LengthFilter(50,500)) // keep protein chains with at least 50 residues
				.filter(new ChainIdFilter<Point3d[]>(chainIdsBc)) // calculate feature vectors for chains in the training set only
		        .mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(1))) // add new chain smoother here ...
//	       	    .mapToPair(new ChainToSequenceFeatureVectorMapper(new AngleSequenceFingerprint()))
//	       	    .mapToPair(new ChainToSequenceFeatureVectorMapper(new DCT1DSequenceFingerprint()))
	       	    .mapToPair(new ChainToSequenceFeatureVectorMapper(new EndToEndDistanceSequenceFingerprint()))
	       	    .cache();
      
        // broadcast feature vectors
        List<Tuple2<String,SequenceFeatureInterface<?>>>  bc =  features.collect(); // return results to master node     
		final Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> featureVectorsBc = sc.broadcast(bc);
        int numVectors = bc.size();
		
		// broadcast list of chain ids that have feature vectors
		List<String> availableChainIds = features
				.keys()
				.collect();
		final Broadcast<List<String>> availableChainIdsBc = sc.broadcast(availableChainIds);
		
		long t2 = System.nanoTime();
		
		// calculate Jaccard Index and join with TM metrics
//	    List<Tuple2<String, Tuple2<Float, String>>> results = pairs
//				.filter(new ChainIdPairFilter(availableChainIdsBc)) // only keep pairs that have feature vectors available
//				.mapToPair(new ChainIdToIndexMapper(availableChainIdsBc)) // map chain ids to indices into feature vector
////				.mapToPair(new LCSFeatureIndexP3(featureVectorsBc,0))
//				.mapToPair(new SmithWatermanP3(featureVectorsBc,0))
//				.mapToPair(new SmithWatermanGotohP3(featureVectorsBc,0))
//				.mapToPair(new JaccardScoreMapperP3(featureVectorsBc))
//				.mapToPair(new LevenshteinMapperP3(featureVectorsBc))
//				.join(trainingData) // join with TM metrics from the input file
//				.sortByKey()
//				.collect();
	    
		sc.stop();
		sc.close();

		
//		int numPairs = results.size();
//		
//		// write results to .csv file
//		PrintWriter writer = new PrintWriter(outputFileName);
//		writeToCsv2(writer, results);
//		writer.close();
//		
//		printStatistics(results);

		long t3 = System.nanoTime();


//     	System.out.println("protein chains     : " + numVectors);
//		System.out.println("total pairs        : " + numPairs);
//		System.out.println();
//		System.out.println("calculate features : " + (t2-t1)/1E9 + " s");
//		System.out.println("compare pairs      : " + (t3-t2)/1E9 + " s");
//		System.out.println("total time         : " + (t3-t1)/1E9 + " s");
//		System.out.println("time per pair      : " + ((t3-t1)/numPairs)  + " ns");
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
	
	/**
	 * 
	 * @param writer
	 * @param list
	 */
	private static void printStatistics(List<Tuple2<String, Tuple2<Float, String>>> joinedResults) {	
		System.out.printf("%7s %7s %7s %7s %7s %7s %7s %7s", "F", "TP", "FN", "TN", "FP", "SENS", "SPEC", "F1");
		System.out.println();
		float tmThreshold = 0.5f;
		for (float f = 0.1f; f < 0.8f; f+= 0.05f) {
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

