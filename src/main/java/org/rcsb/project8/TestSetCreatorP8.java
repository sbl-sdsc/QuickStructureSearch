package org.rcsb.project8;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.hadoop.io.HadoopToSimplePolymerChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralSimilarity.ChainPairToTmMapper;
import org.rcsb.structuralSimilarity.GapFilter;

import scala.Tuple2;

/**
 * This class creates Test Set for future comparing of FatCat TM score
 * The test set will contain different alignment percentage (alignment/shorterChain)
 * 
 * @author  Chris Li, Peter Rose
 */
public class TestSetCreatorP8 { 
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 50;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0]; 
		String outputFilePath = args[1];
		int nPairs = Integer.parseInt(args[2]);
		int seed = Integer.parseInt(args[3]);
		
		long t1 = System.nanoTime();
		TestSetCreatorP8 creator = new TestSetCreatorP8();
		creator.run(sequenceFileName, outputFilePath, nPairs, seed);
		System.out.println("Running Time		: " + ((System.nanoTime()-t1)/1E9) + " s");
	}

	private void run(String path, String outputFilePath, int nPairs, int seed) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName("1" + this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		// Step 1. calculate <pdbId.chainId, feature vector> pairs		
        List<Tuple2<String, Point3d[]>> chains = sc
				.sequenceFile(path, Text.class, SimplePolymerChain.class, NUM_THREADS)  // read protein chains
				.mapToPair(new HadoopToSimplePolymerChainMapper()) // convert input to <pdbId.chainId, protein sequence> pairs
				.filter(t -> t._2.isProtein())
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))	
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				// .filter(new LengthFilter(50,500)) // keep protein chains with 50 - 500 residues
				.collect(); // return results to master node
                
		// Step 2.  broadcast feature vectors to all nodes
		final Broadcast<List<Tuple2<String,Point3d[]>>> chainsBc = sc.broadcast(chains);
		int nChains = chains.size();
		Random r = new Random(seed);
		
		PrintWriter writer1 = new PrintWriter(outputFilePath + "testSet0-20.csv");
		PrintWriter writer2 = new PrintWriter(outputFilePath + "testSet20-40.csv");
		PrintWriter writer3 = new PrintWriter(outputFilePath + "testSet40-60.csv");
		PrintWriter writer4 = new PrintWriter(outputFilePath + "testSet60-80.csv");
		PrintWriter writer5 = new PrintWriter(outputFilePath + "testSet80-100.csv");

						
        // Step 3. map through all pairs for TM score
		for (int i = 0; i < nPairs; i+=BATCH_SIZE) {   
			List<Tuple2<Integer,Integer>> pairs = randomPairs(nChains, BATCH_SIZE, r.nextLong());

			List<Tuple2<String, Float[]>> list = sc
					.parallelizePairs(pairs, NUM_THREADS*NUM_TASKS_PER_THREAD) // distribute data
					//.filter(new ChainPairLengthFilter(chainsBc, 0.5, 1.0)) // restrict the difference in chain length
					.mapToPair(new ChainPairToTmMapper(chainsBc)) // maps pairs of chain id indices to chain id, TM score pairs
					.collect();	// copy result to master node
			// write results to .csv file
			divideSet(writer1,writer2,writer3,writer4,writer5, list);
		}
		
		writer1.close();
		writer2.close();
		writer3.close();
		writer4.close();
		writer5.close();
		
		sc.stop();
		sc.close();

		System.out.println("protein chains     	: " + nChains);
		System.out.println("ramdom pairs       	: " + nPairs);
	}

	/**
	 * Based on the alignment coverage, write to different files
	 * @param writer1
	 * @param writer2
	 * @param writer3
	 * @param writer4
	 * @param writer5
	 * @param list
	 */
	private static void divideSet(PrintWriter writer1, PrintWriter writer2,
			PrintWriter writer3, PrintWriter writer4, PrintWriter writer5,
			List<Tuple2<String, Float[]>> list) {
		for (Tuple2<String, Float[]> t : list) {
			double precent = Math.max(t._2[4],t._2[5]);
			if (precent <= 20)
				writeToCsv(writer1,t);
			else if (precent <= 40)
				writeToCsv(writer2,t);
			else if (precent <= 60)
				writeToCsv(writer3,t);
			else if (precent <= 80)
				writeToCsv(writer4,t);
			else 
				writeToCsv(writer5,t);
		}
	}
	
	/**
	 * write the testset to csv file
	 * @param writer
	 * @param t
	 */
	private static void writeToCsv(PrintWriter writer, Tuple2<String, Float[]> t) {
		writer.print(t._1);
		for (Float f: t._2) {
			writer.print(",");
			writer.print(f);
		}
		writer.println();
		writer.flush();
	}

	/**
	 * Returns random pairs of indices for the pairwise comparison.
	 * @param n
	 * @param nPairs
	 * @param seed
	 * @return
	 */
	private List<Tuple2<Integer, Integer>> randomPairs(int n, int nPairs, long seed) {
		if (nPairs > BATCH_SIZE)
			nPairs = BATCH_SIZE;
		Random r = new Random(seed);
		Set<Tuple2<Integer,Integer>> set = new HashSet<>(nPairs);

		int i = 0;
		int t = 0;
		while (i < nPairs && t < (3 * nPairs)) {
			t++;
			int j = r.nextInt(n);
			int k = r.nextInt(n);
			if (j == k) {
				continue;
			}
			Tuple2<Integer,Integer> tuple = new Tuple2<>(j,k);
			Tuple2<Integer,Integer> tuple2 = new Tuple2<>(k,j);
			if (!set.contains(tuple) && ! set.contains(tuple2)) {
			    set.add(tuple);
			    i++;
			}
		}
		return new ArrayList<Tuple2<Integer,Integer>>(set);
	}
}

