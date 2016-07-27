package org.rcsb.projectec;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.Accumulable;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.project10.BinaryClassificationMapper;
import org.rcsb.project10.WritableSegment;
import org.rcsb.project3.AlignmentAlgorithmInterface;
import org.rcsb.project3.ChainToSequenceFeatureVectorMapper;
import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.project3.SequenceFingerprint;
import org.rcsb.project3.SmithWatermanGotohMapperP3;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/*
 * An attempt at auto-testing many different LibraryFingerprints
 */
public class AutomatedTesting3 {

	private static int NUM_TASKS_PER_THREAD = 3;
	static SparkConf conf = new SparkConf().setMaster("local[*]").setAppName("AutomatedTesting2");
	//
	static JavaSparkContext sc = new JavaSparkContext(conf);

	public static void main(String[] args) throws IOException {
		// Fingerprint Benchmark

		// frag length
		int startLength = 9;
		int endLength = 11;

		double startRmsdThreshold = 3.5;
		double endRmsdThreshold = 5.0;

		String chainsDir = args[0];
		String benchmarkDir = args[1];
		String resultsDir = args[2];

		System.out.println("Chain s       : " + chainsDir);
		System.out.println("Benchmark data: " + benchmarkDir);

		AlignmentAlgorithmInterface algorithm = new SmithWatermanGotohMapperP3();

		JavaRDD<String[]> benchmarkData = sc.textFile(benchmarkDir, sc.defaultParallelism() * NUM_TASKS_PER_THREAD) // read
				// files
				.sample(false, 0.05)
				.map(s -> s.split(","))// split each line into a list items
				// .filter(s -> Double.parseDouble(s[5])>49)
				.cache();

		Set<String> chainIds = new HashSet<>(
				benchmarkData.flatMap(t -> Arrays.asList(t[0], t[1])).distinct().collect());
		JavaPairRDD<String,Point3d[]> coords = sc
				.sequenceFile(chainsDir, Text.class, WritableSegment.class) // read
																			// file
																			// with
																			// chains
				.mapToPair(t -> new Tuple2<String, WritableSegment>(new String(t._1.toString()),
						new WritableSegment(t._2))) // make a copy of the data
				.filter(t -> chainIds.contains(t._1)) // only read chains
														// required for training
														// set
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()));
		
		for (int length = startLength; length <= endLength; length++) {
			for (double rmsdThresh = startRmsdThreshold; rmsdThresh <= endRmsdThreshold; rmsdThresh += .5) {

				// make lib and fingerprint
				List<Point3d[]> library = getLib(chainsDir, length, rmsdThresh);
				double[][] rmsdArray = getRmsdArray(library);
				SequenceFingerprint fingerprint = new LibraryFingerprintMinSim(library, rmsdArray, rmsdThresh);

				AutomatedTesting3 at = new AutomatedTesting3();

				// create unique results directory name
				String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
				resultsDir += File.separatorChar + fingerprint.getName() + "_" + algorithm.getName() + "_" + timeStamp
						+ ".csv";
				System.out.println("Results       : " + resultsDir);
				// run the benchmark
				long start = System.nanoTime();

				at.run(chainsDir, benchmarkDir, resultsDir, fingerprint, algorithm, coords, benchmarkData);
				long end = System.nanoTime();
				double time = (end - start) / 1E9;

				System.out.println("Time          : " + time + " seconds");

				// split input lines of .csv files into chainId1,chainId2, and
				// alignment metrics
				JavaPairRDD<String, Double> tmScores = sc.textFile(benchmarkDir, sc.defaultParallelism()) // read
																											// files
																											// as
																											// lines
						.map(t -> t.split(",")) // split String into String[]
						.mapToPair(t -> new Tuple2<String, Double>(t[0] + "_" + t[1], Double.parseDouble(t[2])));

				// System.out.println("Benchmark size: " + tmScores.count());
				// tmScores.foreach(t -> System.out.println(t));

				// split input lines of .csv files into chainId1,chainId2, and
				// alignment metrics
				JavaPairRDD<String, Double> fpScores = sc.textFile(resultsDir, sc.defaultParallelism()) // read
																		// lines
						.map(t -> t.split(",")) // split String into String[]
						.mapToPair(t -> new Tuple2<String, Double>(t[0] + "_" + t[1], Double.parseDouble(t[2])));

				// System.out.println("Results size: " + fpScores.count());
				// fpScores.foreach(t -> System.out.println(t));

				JavaPairRDD<String, Tuple2<Double, Double>> jScores = tmScores.join(fpScores);

				// calculate binary classification
				double threshold = 0.5;
				JavaPairRDD<String, String> classificationResults = jScores
						.mapToPair(new BinaryClassificationMapper(threshold)).cache();
				//classificationResults.foreach(t -> System.out.println(t));

				// calculate binary classification statistics
				// see: https://en.wikipedia.org/wiki/Precision_and_recall
				long tp = classificationResults.filter(t -> t._2.equals("TP")).count();
				long tn = classificationResults.filter(t -> t._2.equals("TN")).count();
				long fp = classificationResults.filter(t -> t._2.equals("FP")).count();
				long fn = classificationResults.filter(t -> t._2.equals("FN")).count();

				double sensitivity = 1.0 * tp / (double) (tp + fn);
				double specificity = 1.0 * tn / (double) (tn + fp);
				double f1Score = 2.0 * tp / (double) (2 * tp + fp + fn);

				// print binary classification statistics
				System.out.println("Fragment Length: " + length);
				System.out.println("Rmsd Threshold: " + rmsdThresh);
				System.out.println("tp: " + tp);
				System.out.println("tn: " + tn);
				System.out.println("fp: " + fp);
				System.out.println("fn: " + fn);
				System.out.println("Sensitivity: " + sensitivity);
				System.out.println("Specificity: " + specificity);
				System.out.println("F1 score: " + f1Score);

			}
		}

		sc.stop();
		sc.close();
	}

	private static List<Point3d[]> getLib(String chainFile, int fragmentSize, double rmsdThreshold) throws IOException {

		// setup spark
		// SparkConf conf = new
		// SparkConf().setMaster("local[*]").setAppName("ArchLibGenerator").set("spark.serializer",
		// "org.apache.spark.serializer.KryoSerializer");
		//
		// JavaSparkContext sc = new JavaSparkContext(conf);

		// read protein chains and cut into fragments (sliding window approach)
		JavaRDD<Point3d[]> fragments = sc.sequenceFile(chainFile, Text.class, WritableSegment.class)
				.map(t -> t._2.getCoordinates()) // get the coordinates of the
													// protein chains
				.repartition(1) // create a single partition to generate a
								// single fragment library (this cannot be done
								// in parallel!)
				.flatMap(new FlatMapToFragments(fragmentSize)); // flatmap to
																// fragments

		// Create library of fragment archetypes
		List<Point3d[]> prelib = new ArrayList<>();
		List<Point3d[]> lib = Collections.synchronizedList(prelib);
		Accumulable<List<Point3d[]>, Point3d[]> accLib = new Accumulable<>(lib, new AccumuableListPR(rmsdThreshold));

		fragments.foreach(t -> accLib.add(t));

		List<Point3d[]> archetypes = accLib.value();
		return archetypes;
		// System.out.println(archetypes.size());

	}

	public static double[][] getRmsdArray(List<Point3d[]> lib) {
		double[][] rmsdArray = new double[lib.size()][lib.size()];
		SuperPositionQCP qcp = new SuperPositionQCP(true);

		for (int i = 0; i < lib.size(); i++) {
			for (int j = 0; j < lib.size(); j++) {
				qcp.set(lib.get(i), lib.get(j));
				double rmsd = qcp.getRmsd();
				rmsdArray[i][j] = rmsd;
			}
		}

		return rmsdArray;
	}

	private void run(String chainsDir, String benchmarkDir, String resultsDir, SequenceFingerprint fingerprint,
			AlignmentAlgorithmInterface comparisionAlgorithm, JavaPairRDD<String,Point3d[]>  coords, JavaRDD<String[]> benchmarkData) {

		// setup spark
		// SparkConf conf = new
		// SparkConf().setMaster("local[*]").setAppName(this.getClass().getSimpleName());
		//
		// JavaSparkContext sc = new JavaSparkContext(conf);

		// split input lines of .csv files into chainId1,chainId2, and alignment
		// metrics

		// create a list of unique chain ids

		// compute feature vectors
	
		Map<String, SequenceFeatureInterface<?>> featureVectors = 	coords
				.mapToPair(new ChainToSequenceFeatureVectorMapper(fingerprint)) // calculate
																				// sequence
																				// order
																				// dependent
																				// features
																				// (fingerprints)
				.collectAsMap(); // convert JavaRDD to a Java Map

		// broadcast feature vectors to all nodes
		final Broadcast<Map<String, SequenceFeatureInterface<?>>> featureVectorsBc = sc.broadcast(featureVectors);

		// calculate pairwise similarity scores using the training examples
		comparisionAlgorithm.setSequence(featureVectorsBc);

		JavaPairRDD<String, Float> scores = benchmarkData.mapToPair(t -> new Tuple2<String, String>(t[0], t[1])) // chain
																													// Id
																													// pairs
																													// from
																													// training
																													// data
				.mapToPair(comparisionAlgorithm); // calculated pairwise
													// similarity

		// map results to .csv format and save to text file
		System.out.println("almost there");
		scores.filter(t -> t != null)
		.map(t -> new String(t._1 + "," + t._2))
		.saveAsTextFile(resultsDir);
		// System.out.println("aight");
		// a.foreach(f-> System.out.println(f));
		// System.out.println("finally");

		// terminate Spark
	}

}
