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

/**
 * This class print the result of large amount of comparison of origin FatCat and new FatCat
 * This class is used for testing
 * 
 * @author Chris Li
 */
public class Project4ResultPrinter {
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 50;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0]; 
		String outputFileName1 = args[1];
		String outputFileName2 = args[2];
		int nPairs = Integer.parseInt(args[3]);
		int seed = Integer.parseInt(args[4]);
		
		long t1 = System.nanoTime();
		Project4ResultPrinter creator = new Project4ResultPrinter();
		creator.run(sequenceFileName, outputFileName1, outputFileName2, nPairs, seed);
		System.out.println("Overall Running Time	: " + ((System.nanoTime()-t1)/1E9) + " s");
	}

	private void run(String path, String outputFileName1, String outputFileName2, int nPairs, int seed) throws FileNotFoundException {
		Random r = new Random(seed);
		PrintWriter writer1 = new PrintWriter(outputFileName1);
		writer1.print("Protein Length");
		writer1.print(",");
		writer1.print("Test Pairs");
		writer1.print(",");
		writer1.print("Orign Time");
		writer1.print(",");
		writer1.print("Revise Time");
		writer1.print(",");
		writer1.print("Speed Up");
		writer1.print(",");
		writer1.print("Score");
		writer1.println();
		writer1.flush();
		
		PrintWriter writer2 = new PrintWriter(outputFileName2);
		writer2.print("Protein1 Id");
		writer2.print(",");
		writer2.print("Protein2 Id");
		writer2.print(",");
		writer2.print("Origin");
		writer2.print(",");
		writer2.print("Revise");
		writer2.println();
		writer2.flush();
		
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
		while (proteinLength < 2000) {
			if (proteinLength < 500) {
				proteinLength += 50;
				range = 10;
			} else if (proteinLength < 1000) {
				proteinLength += 500;
				range = 100;
			} else {
				proteinLength += 500;
				range = 200;
			}
			List<Tuple2<String, Point3d[]>> chains = listOf(sc, input, proteinLength, range);

			final Broadcast<List<Tuple2<String,Point3d[]>>> chainsBc = sc.broadcast(chains);

			int nChains = chains.size();

			if (!(nChains > 0)) {
				System.out.println("Length of " + proteinLength + " has no protein");
				continue;
			}
			
			List<Tuple2<String, Float[]>> l1 = new ArrayList<Tuple2<String, Float[]>>();
			List<Tuple2<String, Float[]>> l2 = new ArrayList<Tuple2<String, Float[]>>();
			long time1 = 0;
			long time2 = 0;
			
			int PairsLeft = nPairs;
			while (PairsLeft > 0) {

				List<Tuple2<Integer,Integer>> pairs = randomPairs(nChains, PairsLeft, r.nextLong());
				PairsLeft -= BATCH_SIZE;

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
				
				time1 += endTime1 - startTime1;
				time2 += endTime2 - startTime2;
				l1.addAll(list1);
				l2.addAll(list2);
			}

			double score = score(l1,l2);
			
			writeToCsv1(writer1, proteinLength, l1.size(), time1, time2, score);
			
			writeToCsv2(writer2, l1, l2);
		}

		writer1.close();
		writer2.close();
		
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
	private static void writeToCsv1(PrintWriter writer, int proteinLength, int pairs, long time1, long time2, double score) {
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
	
	private static void writeToCsv2(PrintWriter writer, List<Tuple2<String, Float[]>> list1, List<Tuple2<String, Float[]>> list2) {
		for (int i = 0; i < list1.size(); i++) {
			Tuple2<String, Float[]> t1 = list1.get(i);
			Tuple2<String, Float[]> t2 = list2.get(i);
			if (t1._1.equals(t2._1)) {
				writer.print(t1._1);
				for (int j = 0; j < t1._2.length; j++) {
					writer.print(",");
					writer.print(t1._2[j]);
					writer.print(",");
					writer.print(t2._2[j]);
				}
				writer.println();
			}
		}
		writer.flush();
	}
	
	/**
	 * Returns random pairs of indices for the pairwise comparison.
	 * @param n number of feature vectors
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
