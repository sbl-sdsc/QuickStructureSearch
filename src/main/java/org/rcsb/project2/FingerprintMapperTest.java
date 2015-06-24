package org.rcsb.project2;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import javax.vecmath.Point3d;

//
import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.jblas.util.Random;
import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.structuralSimilarity.SeqToChainMapper;

import scala.Tuple2;

public class FingerprintMapperTest {
	private static final int NUM_THREADS = 4;
	private static final int NUM_TASKS_PER_THREAD = 3;

	public static void main(String[] args) {
		String path = args[0];
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]").setAppName(
				FingerprintMapper_KevinWu.class.getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);
		long start = System.nanoTime();
		List<Tuple2<String, Point3d[]>> list = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS * NUM_TASKS_PER_THREAD)
				// .sample(false, 0.002)// sample
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate array> pairs
				// .filter(new GapFilter(0, 0)) // filter chains with zero gap length and zero gaps
				.filter(new org.apache.spark.api.java.function.Function<Tuple2<String, Point3d[]>, Boolean>() {
					@Override
					public Boolean call(Tuple2<String, Point3d[]> v1) throws Exception {
						return v1._1().toUpperCase().equals("1SRR.A") || v1._1().toUpperCase().equals("1T37.A");
					}
				})// filter test
				.collect();
		List<Integer> valid = new ArrayList<>();
		for (int i = 0; i < list.size() - 1; i++)
			if (list.get(i)._2.length == list.get(i + 1)._2.length)
				valid.add(i);

		valid.parallelStream().forEach(new Consumer<Integer>() {
			@Override
			public void accept(Integer i) {
				Point3d[] chain1 = list.get(i)._2;
				Point3d[] chain2 = list.get(i + 1)._2;
				String chainId1 = list.get(i)._1;
				String chainId2 = list.get(i + 1)._1;
				class SeqFeatTest implements SequenceFeatureInterface<Integer> {
					private static final int m = 5;
					private final List<Integer> rand = new ArrayList<>();

					@Override
					public double similarity(SequenceFeatureInterface<Integer> sequence2, int i, int j) {
						return 1 / (1d + (get(i).intValue() ^ sequence2.get(j).intValue()));
					}

					@Override
					public boolean identity(SequenceFeatureInterface<Integer> sequence2, int i, int j) {
						return get(i) == sequence2.get(j);
					}

					@Override
					public Integer[] getSequence() {
						return rand.toArray(new Integer[rand.size()]);
					}

					@Override
					public Integer get(int index) {
						while (rand.size() <= index) {
							rand.add(Random.nextInt(m));
						}
						return rand.get(index);
					}

					@Override
					public int length() {
						return rand.size();
					}

					@Override
					public String toString(int index) {
						return get(i).toString();
					}

					@Override
					public double todouble(int index) {
						// TODO Auto-generated method stub
						return 0;
					}
				}
				StringBuilder test = new StringBuilder(chainId1 + ", " + chainId2);
				test.append(System.lineSeparator());
				test.append(FingerprintMapper_KevinWu.align(chain1, chain2, new SeqFeatTest(), new SeqFeatTest()));
				test.append(System.lineSeparator());
				test.append(System.lineSeparator());
				test.append(FingerprintMapper_KevinWu.align(chain1, chain2));
				System.out.println(test);
			}
		});
		sc.close();
		System.out.println("Time: " + (System.nanoTime() - start) / 1E9 + " sec.");
	}
}
