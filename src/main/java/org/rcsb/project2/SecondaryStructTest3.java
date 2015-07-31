package org.rcsb.project2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.structuralSimilarity.ChainSmootherMapper;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;

import scala.Tuple2;

public class SecondaryStructTest3 {
	private static final int NUM_THREADS = 4;
	private static final int NUM_TASKS_PER_THREAD = 3;
	private static final boolean smooth = true;

	public static void main(String[] args) throws IOException {
		Set<String> needed = new HashSet<>();
		Map<String, SecondaryStruct> pts = new HashMap<>();
		ChainPair[] pairs = new ChainPair[0];
		int N = 0;
		String date;
		date = new SimpleDateFormat("yyyy_MM_dd__hh-mm-ss").format(new Date());
		PrintWriter out = new PrintWriter(new FileWriter("output_" + date + ".csv"));
		try (BufferedReader br = new BufferedReader(new FileReader("data/testsethuge.csv"))) {
			N = Integer.parseInt(br.readLine());
			pairs = new ChainPair[N];
			for (int i = 0; i < N; i++) {
				String[] spl = br.readLine().split(",");
				pairs[i] = new ChainPair();
				pairs[i].setN1(spl[0]);
				pairs[i].setN2(spl[1]);
				pairs[i].setTm(Double.parseDouble(spl[2]));
				System.out.println(pairs[i]);
				needed.add(spl[0]);
				needed.add(spl[1]);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		addAllIntoMap(needed, pts, args[0], smooth);
		for (int i = 0; i < N; i++) {
			System.out.println((i + 1) + " / " + N);
			if (pts.containsKey(pairs[i].getN1()) && pts.containsKey(pairs[i].getN2()))
				out.printf("%s,%s,%.5f,%.5f" + System.lineSeparator(), pairs[i].getN1(), pairs[i].getN2(),
						pairs[i].getTm(),
						SecondaryStructTools.align(pts.get(pairs[i].getN1()), pts.get(pairs[i].getN2())));
		}
		out.close();
	}

	public static void addAllIntoMap(Set<String> needed, Map<String, SecondaryStruct> pts, String path, boolean smooth) {
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(FingerprintMapperTest2.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
		JavaPairRDD<String, Point3d[]> jprdd = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS * NUM_TASKS_PER_THREAD).//
				mapToPair(new HadoopToSimpleChainMapper())// map to hadoop
				.filter(a -> a._2.isProtein())// only proteins
				.filter(a -> a._2.getCoordinates().length > 50)// only length > 50
				.filter(t -> needed.contains(t._1))// names in the set
				.mapToPair(t -> new Tuple2<>(t._1, t._2.getCoordinates()))// string + Point3d[]
				.filter(new GapFilter(0, 0));
		if (smooth)
			jprdd = jprdd.mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(1)));// Smoothing
		pts.putAll(jprdd// s
				.mapToPair(t -> new Tuple2<>(t._1, new SecondaryStruct(t._2, smooth)))// map
				.collectAsMap());
		System.out.println("Needed: " + needed.size());
		System.out.println("Got: " + pts.size());
		sc.close();
	}
}
