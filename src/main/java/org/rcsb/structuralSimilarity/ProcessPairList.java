package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;
import org.rcsb.fingerprints.DCT1DFingerprint;
import org.rcsb.fingerprints.EndToEndDistanceFingerprint;

import scala.Tuple2;

/**
 * 
 * @author  Peter Rose
 */
public class ProcessPairList {  
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 4;
	
	private List<String> chain1 = new ArrayList<>();
	private List<String> chain2 = new ArrayList<>();
	private List<Double> tmScore = new ArrayList<>();

	public static void main(String[] args ) throws IOException
	{
		String oneVsAll = "/Users/peter/Data/ONE_VS_ALL_40/20140919_1105_n_vs_all-40.csv";
		String path = "/Users/peter/Data/PDB_CHAINS/protein_chains_40_20150111_232623.seq";
		String results = "/Users/peter/Data/ProteinSimilarity/FingerPrintAllAgainstAll20150111.csv";

		ProcessPairList aaa = new ProcessPairList();
		aaa.run(path, oneVsAll, results);
	}

	private void run(String path, String oneVsAll, String results) throws IOException {
		readOneVsOne(oneVsAll);
		
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);

		long t1 = System.nanoTime();

		// read file with protein chains and calculate fingerprint vectors
		List<Tuple2<String,Vector>> bc  = sc.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS)  
//				.sample(false, 0.4, 123456)
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate> pairs
				.mapToPair(new ChainToFeatureVectorMapper(new EndToEndDistanceFingerprint())) // calculate fingerprints
	//			.mapToPair(new ChainToFeatureVectorMapper(new DCT1DFingerprint())) // calculate fingerprints
				.collect();


		// broadcast feature vectors to all nodes
		final Broadcast<List<Tuple2<String,Vector>>> vec = sc.broadcast(bc);
		int numVectors = bc.size();
		System.out.println("Vectors: " + numVectors);

		long t2 = System.nanoTime();

		// process pairwise comparisons
		PrintWriter writer = new PrintWriter(results);	

		long numPairs = 0;
		
		List<Tuple2<Integer,Integer>> pairList = getPairList(bc);
		numPairs = pairList.size();
		System.out.println("Pairs: " + numPairs);

		// calculate pairwise scores and filter with score > 0.9
		List<Tuple2<String, Float>> list = sc.parallelizePairs(pairList, NUM_THREADS*NUM_TASKS_PER_THREAD)
				.mapToPair(new PairSimilarityCalculator(vec))
				.collect();	

		// write results to .csv file
		writeToCsv(writer, list);
		long t3 = System.nanoTime();

		writer.close();
		sc.stop();
		sc.close();

		System.out.println("protein chains: " + numVectors);
		System.out.println("pairs         : " + numPairs);
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

	private List<Tuple2<Integer, Integer>> getPairList(List<Tuple2<String,Vector>> bc) {
		Map<String,Integer> map = createChainIdMap(bc);
		
		List<Tuple2<Integer,Integer>> list = new ArrayList<>();
		for (int i = 0; i < this.chain1.size(); i++) {
			Integer index1 = map.get(this.chain1.get(i));
			Integer index2 = map.get(this.chain2.get(i));
			System.out.println(this.chain1.get(i) + "," + this.chain2.get(i));
			System.out.println(index1 + "," + index2);
			if (index1 != null && index2 != null) {
				list.add(new Tuple2<Integer,Integer>(index1, index2));
			}
		}
        return list;
	}
	
	private void readOneVsOne(String fileName) throws IOException {
	      List<String> lines = Files.readAllLines(Paths.get(fileName), Charset.defaultCharset());

	      int count = 0;
	      for (String line : lines) {
	    	  String[] tokens = line.split(",");
	    	  if (count > 0) {
	    		  this.chain1.add(tokens[0]);
	    		  this.chain2.add(tokens[1]);
	    		  this.tmScore.add(Double.parseDouble(tokens[2]));
	    	  }
    		  count++;
	      }
	}
	
	private Map<String,Integer> createChainIdMap(List<Tuple2<String,Vector>> bc) {
		Map<String, Integer> map = new HashMap<String,Integer>();
		for (int i = 0; i < bc.size(); i++) {
			map.put(bc.get(i)._1, i);	
		}
		return map;
	}
}

