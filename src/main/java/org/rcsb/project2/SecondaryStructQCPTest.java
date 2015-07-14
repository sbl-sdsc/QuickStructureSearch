package org.rcsb.project2;

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
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

public class SecondaryStructQCPTest {
	private static final int NUM_THREADS = 4;
	private static final int NUM_TASKS_PER_THREAD = 3;

	public static void main(String[] args) {
		String path = args[0];
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(FingerprintMapperTest2.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
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
		valid.stream().forEach(new Consumer<Integer>() {
			@Override
			public void accept(Integer t) {
				Tuple2<String, SimplePolymerChain> t1 = list.get(t);
				Tuple2<String, SimplePolymerChain> t2 = list.get(t + 1);
				System.out.println(t1._1 + " vs " + t2._1);
				Point3d[] p1 = t1._2.getCoordinates();
				Point3d[] p2 = t2._2.getCoordinates();
				SecondaryStructFinger s1 = new SecondaryStructFinger(p1);
				SecondaryStructFinger s2 = new SecondaryStructFinger(p2);
				if (s1.getAlphaLength() != 0) {
					System.out.println("Alpha Helices");
					System.out.println();
					System.out.println("=(A1, B1)\t=(C1, D1)\t=Segment[E1,F1]");
					System.out.println();
					System.out.println(t1._1);
					s1.printAlphaProjection(0);
					System.out.println(t2._1);
					s2.printAlphaProjection(0);
					System.out.println();
				}
				else {
					System.out.println("no alpha helices");
				}
				if (s1.getBetaLength() != 0) {
					System.out.println("Beta Strands");
					System.out.println(t1._1);
					s1.printBetaProjection(0);
					System.out.println(t2._1);
					s2.printBetaProjection(0);
					System.out.println();
				}
				else {
					System.out.println("no beta strands");
				}
				System.out.println();
			}
		});
		sc.close();
	}

	public static Tuple2<Point3d[], Point3d[]> reduce(Point3d[] p) {
		double[][] dists = SecondaryStructFinger.dists(p);
		Tuple2<int[], int[]> hel = SecondaryStructFinger.alphaHelices(dists, 3);
		int[] hs = hel._1;
		int[] he = hel._2;
		List<Point3d> hp = new ArrayList<>();
		for (int i = 0; i < hs.length; i++) {
			int len = he[i] - hs[i];
			for (int j = 1; j <= len; j++) {
				Point3d po = new Point3d(p[he[i]]);
				po.sub(p[hs[i]]);
				po.scale(j / (double) len);
				po.add(p[hs[i]]);
				hp.add(po);
			}
		}
		Tuple2<int[], int[]> str = SecondaryStructFinger.betaStrands(dists, 4);
		int[] ss = str._1;
		int[] se = str._2;
		List<Point3d> sp = new ArrayList<>();
		for (int i = 0; i < ss.length; i++) {
			int len = se[i] - ss[i];
			for (int j = 1; j <= len; j++) {
				Point3d po = new Point3d(p[se[i]]);
				po.sub(p[ss[i]]);
				po.scale(j / (double) len);
				po.add(p[ss[i]]);
				sp.add(po);
			}
		}
		return new Tuple2<>(hp.toArray(new Point3d[hp.size()]), sp.toArray(new Point3d[sp.size()]));
	}

	public static void compare(Point3d[] x, Point3d[] y) {
		Tuple2<Point3d[], Point3d[]> redx = reduce(x);
		Tuple2<Point3d[], Point3d[]> redy = reduce(y);
		SuperPositionQCP q = new SuperPositionQCP();
		q.set(redx._1, redy._1);
		double a = q.getRmsd();
		q = new SuperPositionQCP();
		q.set(redx._2, redy._2);
		double b = q.getRmsd();
		q = new SuperPositionQCP();
		q.set(x, y);
		System.out.printf("a: %.3f", a);
		System.out.printf("b: %.3f", b);
		System.out.printf("RMSD: %.3f", q.getRmsd());
	}
}
