package org.rcsb.project1;

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
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;

import scala.Tuple2;

/**
 * This class creates structural alignments between random protein chain pairs 
 * using jFatCAT and scores the alignments with the TM score
 * 
 * @author  Peter Rose
 */
public class RandomFragmentCreator { 
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 2; // Spark recommends 2-3 tasks per thread
	private static int BATCH_SIZE = 10000;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0]; 
		String outputFileName = args[1];
		int nPairs = Integer.parseInt(args[2]);
		int seed = Integer.parseInt(args[3]);
		int length = Integer.parseInt(args[4]);
		
		long t1 = System.nanoTime();
		RandomFragmentCreator creator = new RandomFragmentCreator();
		creator.run(sequenceFileName, outputFileName, nPairs, seed, length);
		System.out.println("Time: " + ((System.nanoTime()-t1)/1E9) + " s");
	}

	private void run(String path, String outputFileName, int nPairs, int seed, int length) throws FileNotFoundException {
		// setup spark
		JavaSparkContext sc = getSparkContext();
		
		// Step 1. get <pdbId.chainId, coordinate> pairs
        List<Tuple2<String, Point3d[]>> chains = sc
				.sequenceFile(path, Text.class, SimplePolymerChain.class, NUM_THREADS)  // read protein chains
				.sample(false, 0.1, 123456) // use only a random fraction, i.e., 10%
				.mapToPair(new HadoopToSimplePolymerChainMapper()) // convert input to <pdbId.chainId, SimplePolymerChain> pairs
				.filter(t -> t._2.isProtein())
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(new LengthFilter(50,500)) // keep protein chains with 50 - 500 residues
				.collect(); // return results to master node

		// Step 2.  broadcast coordinates to all nodes
		final Broadcast<List<Tuple2<String,Point3d[]>>> chainsBc = sc.broadcast(chains);
		int nChains = chains.size();

		Random r = new Random(seed);

		PrintWriter writer = new PrintWriter(outputFileName);
		writer.println("PdbId.ChainId1, PdbId.ChainId2,start1,start2,cRMSD,dRMSD,adRMSD,deltacRMSDdRMSD,cRMSDTime,dRMSDTime,adRMSDTime,len1,len2,deltaLen");
		
		for (int i = 0; i < nPairs; i+=BATCH_SIZE) {
			List<Tuple2<Integer,Integer>> pairs = randomPairs(nChains, BATCH_SIZE, r.nextLong());

			List<Tuple2<String, Double[]>> list = sc
					.parallelizePairs(pairs, NUM_THREADS*NUM_TASKS_PER_THREAD) // distribute data
					.mapToPair(new RandomFragmentMapper(chainsBc, length, r.nextInt()))
	//				.filter(t -> t._2[0] < 5.0) // keep only data points with cRMSD < 5 (cRMSD is 0th element in Double[])
					.collect();	// copy result to master node

			// write results to .csv file
			writeToCsv(writer, list);
		}

		writer.close();
		sc.stop();
		sc.close();

		System.out.println("protein chains     : " + nChains);
		System.out.println("ramdom pairs        : " + nPairs);
	}

	
	/**
	 * Writes pairs of chain ids and calculated similarity score to a csv file
	 * @param writer
	 * @param list
	 */
	private static void writeToCsv(PrintWriter writer, List<Tuple2<String, Double[]>> list) {
		for (Tuple2<String, Double[]> t : list) {
			writer.print(t._1);
			for (Double f: t._2) {
				writer.print(",");
			    writer.print(f);
			}
			writer.println();
		}
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
	
	private JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g")
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		return sc;
	}

}

