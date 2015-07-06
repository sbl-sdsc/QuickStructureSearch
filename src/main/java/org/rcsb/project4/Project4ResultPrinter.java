package org.rcsb.project4;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
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
import org.rcsb.structuralSimilarity.ChainPairToTmMapper;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;
import org.rcsb.structuralSimilarity.SeqToChainMapper;

import scala.Tuple2;

public class Project4ResultPrinter {
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 100;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0]; 
		String outputFileName = args[1];
		int nPairs = Integer.parseInt(args[2]);
		int seed = Integer.parseInt(args[3]);
		
		long t1 = System.nanoTime();
		Project4ResultPrinter creator = new Project4ResultPrinter();
		creator.run(sequenceFileName, outputFileName, nPairs, seed);
		System.out.println("Overall Running Time	: " + ((System.nanoTime()-t1)/1E9) + " s");
	}

	private void run(String path, String outputFileName, int nPairs, int seed) throws FileNotFoundException {
		Random r = new Random(seed);
		PrintWriter writer = new PrintWriter(outputFileName);
		writer.print("Protein Length");
		writer.print(",");
		writer.print("Test Pairs");
		writer.print(",");
		writer.print("Orign Time");
		writer.print(",");
		writer.print("Revise Time");
		writer.print(",");
		writer.print("Speed Up");
		writer.print(",");
		writer.print("Score");
		writer.println();
		writer.flush();
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName("1" + this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g");
		
		JavaSparkContext sc = new JavaSparkContext(conf);
		
		JavaPairRDD<String, Point3d[]> input = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS)  // read protein chains
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate> pairs
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.cache();
		
		int proteinLength = 0;
		int range = 10;
		while (proteinLength < 2500) {
			if (proteinLength < 500) {
				proteinLength += 50;
				range = 10;
			} else {
				proteinLength += 500;
				range = 100;
			}
			List<Tuple2<String, Point3d[]>> chains = listOf(sc, input, proteinLength, range);

			final Broadcast<List<Tuple2<String,Point3d[]>>> chainsBc = sc.broadcast(chains);

			int nChains = chains.size();

			if (!(nChains > 0)) {
				System.out.println("Length of " + proteinLength + " has no protein");
				continue;
			}
				
			List<Tuple2<Integer,Integer>> pairs = randomPairs(nChains, BATCH_SIZE, r.nextLong());
			long startTime1 = System.nanoTime();
			List<Tuple2<String, Float[]>> list1 = sc
					.parallelizePairs(pairs, NUM_THREADS*NUM_TASKS_PER_THREAD) // distribute data
					.mapToPair(new ChainPairToTmMapper(chainsBc)) // maps pairs of chain id indices to chain id, TM score pairs
					.collect();	// copy result to master node

			// write results to .csv file
			long endTime1 = System.nanoTime();

			long startTime2 = System.nanoTime();
			List<Tuple2<String, Float[]>> list2 = sc
					.parallelizePairs(pairs, NUM_THREADS*NUM_TASKS_PER_THREAD) // distribute data
					.mapToPair(new ChainPairToTmMapperP4(chainsBc)) // maps pairs of chain id indices to chain id, TM score pairs
					.collect();	// copy result to master node

			// write results to .csv file
			long endTime2 = System.nanoTime();

			double score = score(list1,list2);
			
			writeToCsv(writer, proteinLength, list1.size(), endTime1 - startTime1, endTime2 - startTime2, score);

		}

		writer.close();
		
		sc.stop();
		sc.close();
	}
	
	private double score(List<Tuple2<String, Float[]>> list1, List<Tuple2<String, Float[]>> list2) {
		double score = 0;
		for (int i = 0; i < list1.size(); i ++) {
			float tm1 = list1.get(i)._2[0];
			float tm2 = list2.get(i)._2[0];
			score += Math.abs(tm1 - tm2) * tm1 * 10000;
		}
		return score / list1.size();
	}
	
	private List<Tuple2<String, Point3d[]>> listOf(JavaSparkContext sc, JavaPairRDD<String, Point3d[]> input, int length, int range) {
		List<Tuple2<String, Point3d[]>> chains = input
				.filter(new LengthFilter(length - range, length + range)) // keep protein chains with 50 - 500 residues
				.collect(); // return results to master node
		return chains;
	}
	
	/**
	 * Writes pairs of chain ids and calculated similarity score to a csv file 
	 * @param writer
	 * @param proteinLength
	 * @param score 
	 * @param score2 
	 * @param time 
	 * @param l 
	 */
	private static void writeToCsv(PrintWriter writer, int proteinLength, int pairs, long time1, long time2, double score) {
		writer.print(proteinLength);
		writer.print(",");
		writer.print(pairs);
		writer.print(",");
		writer.print(time1/1E9);
		writer.print(",");
		writer.print(time2/1E9);
		writer.print(",");
		writer.print((time1 - time2)/(time1*1.0) * 100);
		writer.print(",");
		writer.print(score);
		writer.println();
		writer.flush();
	}
	
	/**
	 * Returns random pairs of indices for the pairwise comparison.
	 * @param n number of feature vectors
	 * @return
	 */
	private List<Tuple2<Integer, Integer>> randomPairs(int n, int nPairs, long seed) {
		Random r = new Random(seed);
		Set<Tuple2<Integer,Integer>> set = new HashSet<>(nPairs);

		for (int i = 0; i < nPairs; i++) {
			int j = r.nextInt(n);
			int k = r.nextInt(n);
			if (j == k) {
				continue;
			}

			Tuple2<Integer,Integer> tuple = new Tuple2<>(j,k);
			if (! set.contains(tuple)) {
			    set.add(tuple);
			}
		}
		return new ArrayList<Tuple2<Integer,Integer>>(set);
	}
}
