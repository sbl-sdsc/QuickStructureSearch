package org.rcsb.projectva;
import org.rcsb.project3.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.project10.WritableSegment;
import org.rcsb.project3.AlignmentAlgorithmInterface;
import org.rcsb.project3.ChainToSequenceFeatureVectorMapper;
import org.rcsb.project3.EndToEndDistanceSequenceFingerprint;
//import org.rcsb.project3.MeetMinIndexMapper;
import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.project3.SequenceFingerprint;
//import org.rcsb.projectec.ArchLibGeneratorPR;
//import org.rcsb.projectec.LibraryFingerprint;

import scala.Tuple2;

/**
 * This class reads a benchmark set of protein chain pairs and calculates
 * similarity scores using fingerprinting methods.
 * The scores are saved as .csv files.
 * 
 * @author Peter Rose
 */
public class FingerprintBenchmarkva implements Serializable {
	private static final long serialVersionUID = -8293414734009053770L;	
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args) throws FileNotFoundException {
		if (args.length != 3) {
			System.out.println("Usage: FingerPrintTester.jar chainsDir trainingDir resultsDir");
			System.out.println("  chainsDir: directory with Hadoop sequence files containing protein chains");
			System.out.println("  benchmarkDir: directory with benchmark data");
			System.out.println("  resultsDir: directory with calculated scores");
		}
		
		String chainsDir = args[0];
		String benchmarkDir = args[1];
		String resultsDir = args[2];
		
		System.out.println("Chain s       : " + chainsDir);
		System.out.println("Benchmark data: " + benchmarkDir);
		
		List<Point3d[]> library = ArchLibGeneratorPR.readLibraryFromFile(args[3]);
		// setup fingerprint algorithm
		
//		SequenceFingerprint fingerprint = new EndToEndDistanceSequenceFingerprint();
//		SequenceFingerprint fingerprint = new DCT1DSequenceFingerprint();
		SequenceFingerprint fingerprint = new LibraryFingerprint(library,2.0);
		
		// setup similarity algorithm
//		AlignmentAlgorithmInterface algorithm = new NormalizedCompressionDistanceMapper();
//		AlignmentAlgorithmInterface algorithm = new LevenshteinMapperP3();
	    AlignmentAlgorithmInterface algorithm = new SmithWatermanGotohMapperP3();

		FingerprintBenchmarkva benchmark = new FingerprintBenchmarkva();
		
		// create unique results directory name
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		resultsDir += File.separatorChar + fingerprint.getName() + "_" + algorithm.getName() + "_"+ timeStamp + ".csv";
		System.out.println("Results       : " + resultsDir);
	
		// run the benchmark
		long start = System.nanoTime();
		
		benchmark.run(chainsDir, benchmarkDir, resultsDir, fingerprint, algorithm);
			    
	    long end = System.nanoTime();   
	    double time = (end-start)/1E9;
	    
	    System.out.println("Time          : " + time + " seconds"); 
	}
//
	/**
	 * Runs protein chain similarity calculations for a benchmark set.
	 * @param chainsDir directory with Hadoop sequence files containing protein chains
	 * @param benchmarkDir directory with benchmark data
	 * @param resultsDir directory with calculated score
	 * @param fingerprint fingerprinting method
	 * @param comparisionAlgorithm algorithm for sequence feature vector comparision
	 */
	private void run(String chainsDir, String benchmarkDir, String resultsDir, SequenceFingerprint fingerprint, 
			AlignmentAlgorithmInterface comparisionAlgorithm) {
		
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[*]")
				.setAppName(this.getClass().getSimpleName());	
		
		JavaSparkContext sc = new JavaSparkContext(conf);

		// split input lines of .csv files into chainId1,chainId2, and alignment metrics
		JavaRDD<String[]> benchmarkData = sc
				.textFile(benchmarkDir, sc.defaultParallelism() * NUM_TASKS_PER_THREAD) // read files
				.map(s -> s.split(",")) // split each line into a list items
				.cache();	

		// create a list of unique chain ids
		Set<String> chainIds = new HashSet<>
		      (benchmarkData.flatMap(t -> Arrays.asList(t[0], t[1])).distinct().collect());
		
		// compute feature vectors	    
	    Map<String, SequenceFeatureInterface<?>> featureVectors = sc
	    		.sequenceFile(chainsDir, Text.class, WritableSegment.class) // read file with chains
	    		.mapToPair(t -> new Tuple2<String, WritableSegment> (new String(t._1.toString()), new WritableSegment(t._2)) ) // make a copy of the data
	    		.filter(t -> chainIds.contains(t._1)) // only read chains required for training set
	    		.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates())) // map to chain Id, coordinate pairs
				.mapToPair(new ChainToSequenceFeatureVectorMapper(fingerprint)) // calculate sequence order dependent features (fingerprints)
				.collectAsMap(); // convert JavaRDD to a Java Map

	    // broadcast feature vectors to all nodes
	    final Broadcast<Map<String, SequenceFeatureInterface<?>>> featureVectorsBc = sc
				.broadcast(featureVectors);
	    
	    // calculate pairwise similarity scores using the training examples
	    comparisionAlgorithm.setSequence(featureVectorsBc);

	    JavaPairRDD<String, Float> scores = benchmarkData
	    		.mapToPair(t -> new Tuple2<String, String>(t[0], t[1])) // chain Id pairs from training data
	    		.mapToPair(comparisionAlgorithm); // calculated pairwise similarity
		
		// map results to .csv format and save to text file
		scores
				.filter(t -> t!= null)
		        .map(t -> new String(t._1 + "," + t._2))
		        .saveAsTextFile(resultsDir);

	    // terminate Spark
	    sc.stop();
	    sc.close();
	}
}
