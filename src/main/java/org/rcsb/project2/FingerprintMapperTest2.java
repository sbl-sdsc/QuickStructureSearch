package org.rcsb.project2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.project3.SequenceFeatureInterface;

import scala.Tuple2;

public class FingerprintMapperTest2 {
	private static final int NUM_THREADS = 4;
	private static final int NUM_TASKS_PER_THREAD = 3;

	public static void main(String[] args) throws IOException {
		String[] id1, id2;
		BufferedReader br = new BufferedReader(new FileReader("data/pairs.txt"));
		int N = Integer.parseInt(br.readLine());
		id1 = new String[N];
		id2 = new String[N];
		for (int i = 0; i < N; i++) {
			String[] spl = br.readLine().split(",");
			id1[i] = spl[0];
			id2[i] = spl[1];
			try {
				System.out.println(id1[i]);
				SecondaryStructTools.write(SecondaryStructTools.pull(id1[i]), id1[i]);
				System.out.println(id2[i]);
				SecondaryStructTools.write(SecondaryStructTools.pull(id2[i]), id2[i]);
			}
			catch (Exception e) {
				System.err.println("BAD " + id1[i] + " : " + id2[i]);
				e.printStackTrace();
			}
		}
		br.close();
	}

	public static void run2(String[] args) {
		String path = args[0];
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(FingerprintMapperTest2.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
		List<Tuple2<String, SimplePolymerChain>> list = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS * NUM_TASKS_PER_THREAD)
				.sample(false, 0.02, 2346)// sample
				.mapToPair(new HadoopToSimpleChainMapper())
				// .filter(new GapFilter(0, 0)) // filter chains with zero gap length and zero gaps
				.filter(t -> t._2.isProtein())// filter test
				.collect();
		List<Integer> valid = new ArrayList<>();
		for (int i = 0; i < list.size() - 1; i++)
			if (list.get(i)._2.getCoordinates().length == list.get(i + 1)._2.getCoordinates().length)
				valid.add(i);
		for (int i : valid) {
			System.out.println(list.get(i)._1 + ", " + list.get(i + 1)._1);
		}
		sc.close();
	}

	public static void run1(String[] args) {
		String path = args[0];
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(FingerprintMapperTest2.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
		long start = System.nanoTime();
		List<Tuple2<String, SimplePolymerChain>> list = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS * NUM_TASKS_PER_THREAD)
				.sample(false, 0.002, 2345)// sample
				.mapToPair(new HadoopToSimpleChainMapper())
				// .filter(new GapFilter(0, 0)) // filter chains with zero gap length and zero gaps
				.filter(t -> t._2.isProtein())// filter test
				.collect();
		List<Integer> valid = new ArrayList<>();
		for (int i = 0; i < list.size() - 1; i++)
			if (list.get(i)._2.getCoordinates().length == list.get(i + 1)._2.getCoordinates().length)
				valid.add(i);
		System.out.println("valid: " + valid.size());
		valid.stream().forEach(new Consumer<Integer>() {
			@Override
			public void accept(Integer i) {
				SimplePolymerChain chain1 = list.get(i)._2;
				SimplePolymerChain chain2 = list.get(i + 1)._2;
				String chainId1 = list.get(i)._1;
				String chainId2 = list.get(i + 1)._1;
				class SeqFeatTest implements SequenceFeatureInterface<Integer> {
					private static final int m = 5;
					private final List<Integer> rand = new ArrayList<>();

					@Override
					public double similarity(SequenceFeatureInterface<Integer> sequence2, int i, int j) {
						return 0.04 + 1 / (1d + (get(i).intValue() ^ sequence2.get(j).intValue()));
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
							rand.add((int) (Math.random() * m));
						}
						return rand.get(index);
					}

					@Override
					public int length() {
						return rand.size();
					}

					@Override
					public String toString(int index) {
						return get(index).toString();
					}

					@Override
					public double todouble(int index) {
						return 0;
					}

					@Override
					public Point3d[] getCoords() {
						return null;
					}
				}
				StringBuilder test = new StringBuilder(chainId1 + ", " + chainId2);
				test.append(System.lineSeparator());
				test.append(FingerprintMapper.align(chain1, chain2, new SeqFeatTest(), new SeqFeatTest(), true));
				test.append(System.lineSeparator());
				// test.append(System.lineSeparator());
				// test.append(FingerprintMapper_KevinWu
				// .align(chain1, chain2, new SeqFeatTest(), new SeqFeatTest(), false));
				// test.append(System.lineSeparator());
				test.append(System.lineSeparator());
				System.out.println(test);
			}
		});
		sc.close();
		System.out.println("Time: " + (System.nanoTime() - start) / 1E9 + " sec.");
	}
}
