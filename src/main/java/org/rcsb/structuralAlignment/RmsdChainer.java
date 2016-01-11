package org.rcsb.structuralAlignment;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.regression.LabeledPoint;
import org.rcsb.hadoop.io.HadoopToSimplePolymerChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;

import scala.Tuple2;

/**
 * This class creates structural alignments between random protein chain pairs 
 * using a new algorithm based on the alignment for chain segments
 * 
 * @author  Peter Rose
 */
public class RmsdChainer implements Serializable { 
	private static final long serialVersionUID = 2779213801755875110L;
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0]; 
		int nPairs = Integer.parseInt(args[1]);
		int seed = Integer.parseInt(args[2]);
		
		long t1 = System.nanoTime();
		RmsdChainer creator = new RmsdChainer();
		creator.run(sequenceFileName, nPairs, seed);
		System.out.println("Time: " + ((System.nanoTime()-t1)/1E9) + " s");
	}

	private void run(String path, int nPairs, int seed) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g")
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
				.registerKryoClasses(new Class[]{HadoopToSimplePolymerChainMapper.class});

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		
		// read sequence file and retrieve tuples of PDBId.ChainId, CA Coordinate arrays
		List<Tuple2<String, Point3d[]>> chains = sc
				.sequenceFile(path, Text.class, SimplePolymerChain.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.1, 123456) // use only 10% of the data
				.mapToPair(new HadoopToSimplePolymerChainMapper())
				.filter(t -> t._2.isProtein())
				.map(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
				.filter(new GapFilter(0, 0)) // filter out chains with gaps
				.filter(new LengthFilter(50,500)) // keep protein chains with at least 50 residues
		//		.mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(1))) // add new chain smoother here ...
				.collect();

		// Step 2.  broadcast chains to all nodes
		final Broadcast<List<Tuple2<String,Point3d[]>>> chainsBc = sc.broadcast(chains);
		int nChains = chains.size();

		Random r = new Random(seed);
		
		List<Tuple2<Integer, Integer>> pairs = randomPairs(nChains, nPairs, r.nextLong());

		JavaPairRDD<String, Float> data = sc
				.parallelizePairs(pairs, NUM_THREADS * NUM_TASKS_PER_THREAD)
		//						.mapToPair(new AlignmentMapper(chainsBc))
				.mapToPair(new AlignmentMapper2(chainsBc))
				.filter(s -> s._2 > 0.5) // keep results with a TM score > 0.5
				.cache();

		System.out.println("***** results ******");
		List<Tuple2<String, Float>> hits = data.collect();
		for (Tuple2<String, Float> t: hits) {
			System.out.println(t);
		}
		System.out.println("# hits: " + hits.size());
		
		sc.stop();
		sc.close();

		System.out.println("protein chains     : " + nChains);
		System.out.println("ramdom pairs        : " + nPairs);

	}

	/**
	 * Writes pairs to a csv file
	 * @param list
	 * @param list
	 */
	private static void writeToCsv(PrintWriter writer, List<LabeledPoint> list) {
		for (LabeledPoint p : list) {
			writer.print(p.label());
			writer.print(",");
			writer.println(p.features().toArray()[0]);
		}
		// writer.flush();
	}

	/**
	 * Returns random pairs of indices for the pairwise comparison.
	 * @param n number of feature vectors
	 * @return
	 */
	private List<Tuple2<Integer, Integer>> randomPairs(int n, int nPairs, long seed) {
		Random r = new Random(seed);
		List<Tuple2<Integer,Integer>> list = new ArrayList<>(nPairs);

		for (int i = 0; i < nPairs; i++) {
			int j = r.nextInt(n);
			int k = r.nextInt(n);
			if (j == k) {
				continue;
			}

			Tuple2<Integer,Integer> tuple = new Tuple2<>(j,k);
			if (! list.contains(tuple)) {
			    list.add(tuple);
			}
		}
		return list;
	}
}

