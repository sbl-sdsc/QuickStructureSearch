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
import org.rcsb.project3.AngleSequenceFingerprint;
import org.rcsb.project3.ChainToSequenceFeatureVectorMapper;
import org.rcsb.project3.DCT1DSequenceFingerprint;
import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.project3.SmithWatermanGotohP3;
import org.rcsb.structuralSimilarity.ChainPairToTmMapper;
import org.rcsb.structuralSimilarity.ChainSmootherMapper;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;

import scala.Tuple2;

public class OneAgainstAllP8 {
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 50;
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		if (args.length < 2) {
			System.out.println("Usage: OneAgainstAll.jar proteinId [sequenceFile] outputPath");
			System.out.println("  proteinId: proteinId that will be used to compare with all other proteins");
			System.out.println("  sequenceFile: Hadoop sequence file with protein chain information");
			System.out.println("  outputPath: Path for results from calculation in .csv format.");
		}
		String sequenceFile = "/Users/Chris/Documents/RCSB/Data/Protein_chains/protein_chains_All_20150629_002251.seq";
		String outputPath = args[1];

		if (args.length == 3) {
			sequenceFile = args[1]; 
			outputPath = args[2];
		}
		OneAgainstAllP8 o = new OneAgainstAllP8();
		long st = System.nanoTime();
		o.run(args[0], sequenceFile, outputPath);
		System.out.println("Total Running Time: " + (System.nanoTime()-st)/1E9 + " s");
	}
	
	private void run(String proteinId, String sequenceFile, String outputPath) throws FileNotFoundException {
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName("1" + this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		JavaPairRDD<String, Point3d[]> proteinChains = sc
				.sequenceFile(sequenceFile, Text.class, ArrayWritable.class, NUM_THREADS*NUM_TASKS_PER_THREAD)  // read protein chains
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
		
		// calculate <chainId, feature vector> pairs
        JavaPairRDD<String, SequenceFeatureInterface<?>> features = proteinChains
		        .mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(1))) // add new chain smoother here ...
//	       	    .mapToPair(new ChainToSequenceFeatureVectorMapper(new AngleSequenceFingerprint()))
	       	    .mapToPair(new ChainToSequenceFeatureVectorMapper(new DCT1DSequenceFingerprint()))
//	       	    .mapToPair(new ChainToSequenceFeatureVectorMapper(new EndToEndDistanceSequenceFingerprint()))
	       	    .cache();
        
        // broadcast feature vectors
        List<Tuple2<String,SequenceFeatureInterface<?>>> featureVectors =  features.collect(); // return results to master node     
		
        Integer targetChain = null;
		for (int i = 0; i < featureVectors.size(); i++) {
			Tuple2<String, SequenceFeatureInterface<?>> t = featureVectors.get(i);
			if (t._1.equals(proteinId)) {
				targetChain = i;
				break;
			}
		}
		
		if (targetChain == null) {
			sc.stop();
			sc.close();
			System.out.println("ProteinId " + proteinId + " is not contained in the sequence file.");
			return;
		}
        
		PrintWriter writer = new PrintWriter(outputPath + "OneAgainstAll_" + proteinId + ".csv");
		writer.print("Target Protein");
		writer.print(",");
		writer.print("Protein2");
		writer.print(",");
		writer.println("Score");
		writer.flush();
		
		PrintWriter writer2 = new PrintWriter(outputPath + "OneAgainstAll_" + proteinId + "_TM.csv");
		writer.print("Target Protein");
		writer.print(",");
		writer.print("Protein2");
		writer.print(",");
		writer.println("TM Score");
		writer.flush();

		
		final Integer targetId = targetChain;
        final Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> featureVectorsBc = sc.broadcast(featureVectors);
		
		List<Integer> ChainIndex = new ArrayList<Integer>();
		for (int i = 0; i < featureVectors.size(); i++) {
			ChainIndex.add(i);
		}
		// calculate Jaccard Index and join with TM metrics
	    List<Tuple2<String, Float>> results = sc
	    		.parallelize(ChainIndex)
				.mapToPair(i -> new Tuple2<Integer,Integer>(targetId, i))
//				.mapToPair(new LCSFeatureIndexP3(featureVectorsBc,0))
//				.mapToPair(new SmithWatermanP3(featureVectorsBc,0))
				.mapToPair(new SmithWatermanGotohP3(featureVectorsBc,0))
//				.mapToPair(new JaccardScoreMapperP3(featureVectorsBc))
//				.mapToPair(new LevenshteinMapperP3(featureVectorsBc))
				.sortByKey()
				.collect();
		
	    List<Integer> pass = new ArrayList<Integer>();
	    for (Tuple2<String, Float> t: results) {
	    	if (t._2 < 0.4) {
	    		String[] pros = t._1.split(",");
	    		pass.add(ChainIds.indexOf(pros[1]));
	    	}
	    }
	    targetChain = ChainIds.indexOf(proteinId);
	    
        final Broadcast<List<Tuple2<String,Point3d[]>>> sequence = sc.broadcast(Chains);
		        
        boolean flag = true;
        while(flag) {  
			List<Tuple2<Integer,Integer>> pairs = getPairs(pass, targetChain, BATCH_SIZE);
			System.out.println(pass.size());
			if (pairs.size() < BATCH_SIZE) {
				flag = false;
				if (pairs.size() == 0)
					break;
			}
			
			List<Tuple2<String, Float[]>> list = sc
					.parallelizePairs(pairs, NUM_THREADS*NUM_TASKS_PER_THREAD) // distribute data
					//.filter(new ChainPairLengthFilter(chainsBc, 0.5, 1.0)) // restrict the difference in chain length
					.mapToPair(new ChainPairToTmMapper(sequence)) // maps pairs of chain id indices to chain id, TM score pairs
					.collect();	// copy result to master node
			
			writeToCsv2(writer2,list);
			// write results to .csv file
		}
	    
	    writeToCsv(writer,results);
	    
	    writer.close();
		sc.stop();
		sc.close();
	}

	private void writeToCsv2(PrintWriter writer, List<Tuple2<String, Float[]>> list) {
		for (Tuple2<String, Float[]> t : list) {
			writer.print(t._1);
			for (float score : t._2) {
				writer.print(",");
				writer.print(score);
			}
			writer.println();
		}
		writer.flush();
	}

	private List<Tuple2<Integer, Integer>> getPairs(List<Integer> pass, Integer targetChain, int BATCH_SIZE) {
		List<Tuple2<Integer, Integer>> pairs = new ArrayList<Tuple2<Integer, Integer>>();
		for (int i = 0; i < BATCH_SIZE; i++){
			pairs.add(new Tuple2<Integer, Integer>(targetChain, pass.remove(0)));
			if (pass.size() == 0)
				return pairs;
		}
		return pairs;
	}

	private void writeToCsv(PrintWriter writer, List<Tuple2<String, Float>> results) {
		for (Tuple2<String, Float> t : results) {
			writer.print(t._1);
			writer.print(",");
			writer.printf("%.3f",t._2);
			writer.println();
		}
		writer.flush();
	}
}
