package org.rcsb.project8;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.project3.AlignmentAlgorithmInterface;
import org.rcsb.project3.ChainToSequenceFeatureVectorMapper;
import org.rcsb.project3.EndToEndDistanceDoubleSequenceFingerprint;
import org.rcsb.project3.EndToEndDistanceSequenceFingerprint;
import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.project3.SequenceFingerprint;
import org.rcsb.project3.SmithWatermanGotohP3;
import org.rcsb.project3.SmithWatermanWithGeoComp;
import org.rcsb.project4.ChainPairToTmMapperP4;
import org.rcsb.structuralSimilarity.ChainPairToTmMapper;
import org.rcsb.structuralSimilarity.ChainSmootherMapper;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;

import scala.Tuple2;

/**
 * This class is used of one-against-all computation. It will first run a fingerprint alignment as a threshold before 
 * actual TM score computation.
 * 
 * @author Chris Li
 */
public class OneAgainstAllWithThreshold {
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 50;
	private static double fingerPrintFilter = 0.75;
	// Other fingerPrint: AngleSequenceFingerprint() || DCT1DSequenceFingerprint()
	private static SequenceFingerprint fingerPrint = new EndToEndDistanceDoubleSequenceFingerprint();
	// Other alignment: LCSFeatureIndexP3() || SmithWatermanP3() || JaccardScoreMapperP3() || LevenshteinMapperP3()
	private static AlignmentAlgorithmInterface alignmentAlgorithm = new SmithWatermanWithGeoComp();
	//private static final double[] filters = {0.35,0.5,0.7,0.75};
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		if (args.length < 3) {
			System.out.println("Usage: OneAgainstAll.jar proteinId [sequenceFile] outputPath");
			System.out.println("  proteinId: proteinId that will be used to compare with all other proteins");
			System.out.println("  TM filter: TM score that protein pairs have TM score above will be shown");
			System.out.println("  sequenceFile: Hadoop sequence file with protein chain information");
			System.out.println("  outputPath: Path for results from calculation in .csv format.");
		}
		String sequenceFile = "/Users/Chris/Documents/RCSB/Data/Protein_chains/protein_chains_All_20150629_002251.seq";
		String outputPath = args[1];
		String targetId = "";
		double tmFilter = 0.5;
		if (args.length == 4) {
			targetId = args[0];
			tmFilter = Double.parseDouble(args[1]);
			sequenceFile = args[2]; 
			outputPath = args[3];
		}
		
		if (tmFilter > 0.85)
			fingerPrintFilter = 0.75;
		else if (tmFilter > 0.8)
			fingerPrintFilter = 0.7;
		else if (tmFilter > 0.75)
			fingerPrintFilter = 0.5;
		else if (tmFilter > 0.7)
			fingerPrintFilter = 0.35;
		else fingerPrintFilter = 0.2;
		
		OneAgainstAllWithThreshold o = new OneAgainstAllWithThreshold();
		long st = System.nanoTime();
		o.run(targetId, sequenceFile, outputPath);
		System.out.println("Total Running Time: " + (System.nanoTime()-st)/1E9 + " s");
	}
	
	private void run(String targetProteinId, String sequenceFile, String outputPath) throws FileNotFoundException {
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName("1" + this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		JavaPairRDD<String, Point3d[]> proteinChains = sc
				.sequenceFile(sequenceFile, Text.class, ArrayWritable.class, NUM_THREADS * NUM_TASKS_PER_THREAD)  // read protein chains
				.mapToPair(new HadoopToSimpleChainMapper()) // convert input to <pdbId.chainId, protein sequence> pairs
				.filter(t -> t._2.isProtein())
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))	
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(new LengthFilter(20,3000)) // keep protein chains with 50 - 500 residues
				.cache(); // return results to master node		
		
		List<Tuple2<String, Point3d[]>> Chains = proteinChains.collect();
		        
		List<String> ChainIds = proteinChains
				.keys()
				.collect();
		
		final Integer targetChain = ChainIds.indexOf(targetProteinId);
		
		// check if the target protein is contained in the sequence file
		if (targetChain == -1) {
			sc.stop();
			sc.close();
			System.out.println("ProteinId " + targetProteinId + " is not contained in the sequence file.");
			return;
		}
		
		// 1st step:
		//		calculate <chainId, feature vector> pairs for fingerPrint sequence
        JavaPairRDD<String, SequenceFeatureInterface<?>> features = proteinChains
		        .mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(1))) // add new chain smoother here ...
	       	    .mapToPair(new ChainToSequenceFeatureVectorMapper(fingerPrint))
	       	    .cache();
        
        // broadcast feature vectors
        List<Tuple2<String,SequenceFeatureInterface<?>>> featureVectors =  features.collect(); // return results to master node     
        final Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> featureVectorsBc = sc.broadcast(featureVectors);

		alignmentAlgorithm.setSequence(featureVectorsBc);
		
        final Broadcast<List<Tuple2<String,Point3d[]>>> sequence = sc.broadcast(Chains);

		alignmentAlgorithm.setCoords(sequence);
		
		// ready for pair parallelize
		List<Integer> ChainIndex = new ArrayList<Integer>();
		for (int i = 0; i < featureVectors.size(); i++) {
			ChainIndex.add(i);
		}

		// 2nd step:
		//		calculate alignment score
	    List<Tuple2<String, Float>> results = sc
	    		.parallelize(ChainIndex)
				.mapToPair(i -> new Tuple2<Integer,Integer>(targetChain, i))
				.mapToPair(alignmentAlgorithm)
				.sortByKey()
				.collect();
		
	    // 3rd step:
	    //		use fingerprint filter to get reduce TM computation
	    List<Integer> pass = new ArrayList<Integer>();
	    for (Tuple2<String, Float> t: results) {
	    	if (t._2 >= fingerPrintFilter) {
	    		String[] pros = t._1.split(",");
	    		pass.add(ChainIds.indexOf(pros[1]));
	    	}
	    }
	    
        // output fingerPrint result
		PrintWriter writer = new PrintWriter(outputPath + "OneAgainstAll_" + targetProteinId + "_FingerPrint.csv");
		writer.print("Target Protein");
		writer.print(",");
		writer.print("Protein2");
		writer.print(",");
		writer.println("Score");
		writer.flush();
	    writeFingerPrintResultToCsv(writer,results);
	    writer.close();
		
		PrintWriter writer2 = new PrintWriter(outputPath + "OneAgainstAll_" + targetProteinId + "_TM.csv");
		writer2.print("Target Protein");
		writer2.print(",");
		writer2.print("Protein2");
		writer2.print(",");
		writer2.println("TM Score");
		writer2.flush();
		
        // 4th step:
        //		calculate TM score using FatCat
        boolean flag = true;
        while(flag) { 
        	// run part of the pairs at each time
			List<Tuple2<Integer,Integer>> pairs = generatePairs(pass, targetChain, BATCH_SIZE);
			System.out.println(pass.size());
			if (pairs.size() < BATCH_SIZE) {
				flag = false;
				if (pairs.size() == 0)
					break;
			}
			
			List<Tuple2<String, Float[]>> tmResults = sc
					.parallelizePairs(pairs, NUM_THREADS*NUM_TASKS_PER_THREAD) // distribute data
					.mapToPair(new ChainPairToTmMapperP4(sequence)) // maps pairs of chain id indices to chain id, TM score pairs
					.collect();	// copy result to master node
			
			writeTmScoreToCsv(writer2,tmResults);
        }
	    
	    writer2.close();
	    
		sc.stop();
		sc.close();
	}

	/**
	 * Generate pairs of protein index for future one-against-all computation
	 * @param proteinIds
	 * @param targetProteinId
	 * @param BATCH_SIZE
	 * @return
	 */
	private List<Tuple2<Integer, Integer>> generatePairs(List<Integer> proteinIds, Integer targetProteinId, int BATCH_SIZE) {
		List<Tuple2<Integer, Integer>> pairs = new ArrayList<Tuple2<Integer, Integer>>();
		for (int i = 0; i < BATCH_SIZE; i++){
			if (proteinIds.size() == 0)
				return pairs;
			pairs.add(new Tuple2<Integer, Integer>(targetProteinId, proteinIds.remove(0)));
		}
		return pairs;
	}

	/**
	 * Write the result of fingerprint alignment to file
	 * @param writer
	 * @param results
	 */
	private void writeFingerPrintResultToCsv(PrintWriter writer, List<Tuple2<String, Float>> results) {
		for (Tuple2<String, Float> t : results) {
			writer.print(t._1);
			writer.print(",");
			writer.printf("%.3f",t._2);
			writer.println();
		}
		writer.flush();
	}
	
	/**
	 * Write the result of TM score to file
	 * @param writer
	 * @param results
	 */
	private void writeTmScoreToCsv(PrintWriter writer, List<Tuple2<String, Float[]>> results) {
		for (Tuple2<String, Float[]> t : results) {
			writer.print(t._1);
			for (float score : t._2) {
				writer.print(",");
				writer.print(score);
			}
			writer.println();
		}
		writer.flush();
	}
}
