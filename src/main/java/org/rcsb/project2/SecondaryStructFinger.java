package org.rcsb.project2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import javax.vecmath.Point3d;

public class SecondaryStructFinger {
	private static final int NUM = 4;
	// private static final int NUM_THREADS = 4;
	// private static final int NUM_TASKS_PER_THREAD = 3;

	private double[][] dists;

	public SecondaryStructFinger(Point3d[] pts) {
		dists = dists(pts);
	}

	public double[] get(int i) {
		return dists[i];
	}

	public double[][] getRange(int i, int j) {
		if (i < 0 || i > dists.length || j < 0 || j > dists.length || i > j)
			return null;
		return Arrays.copyOfRange(dists, i, j);
	}

	public void printHelices() {
		int[] hel = alphaHelices(dists, 0);
		if (hel.length % 2 != 0)
			new Exception("Helices length is not even").printStackTrace();
		for (int i = 0; i < hel.length; i += 2) {
			System.out.printf("%d-%d" + System.lineSeparator(), hel[i] - 2, hel[i + 1]);
		}
	}

	public static void main(String[] args) {
		// String path = args[0];
		// SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
		// .setAppName(FingerprintMapperTest2.class.getSimpleName())
		// .set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		// JavaSparkContext sc = new JavaSparkContext(conf);
		// List<Tuple2<String, SimplePolymerChain>> list = sc
		// .sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS * NUM_TASKS_PER_THREAD)
		// // .sample(false, 0.002, 2345)// sample
		// .mapToPair(new HadoopToSimpleChainMapper())
		// // .filter(new GapFilter(0, 0)) // filter chains with zero gap length and zero gaps
		// .filter(t -> t._1.equals("3WST.G"))
		// // .filter(t -> t._2.isProtein())// filter test
		// .collect();
		// long t = System.nanoTime();
		// list.parallelStream().forEach(
		// a -> System.out.println(a._1 + System.lineSeparator() + distsToString(dists(a._2.getCoordinates()))));
		// System.out.println((System.nanoTime() - t) / 1E9);
		// Point3d[] pts = list.get(0)._2.getCoordinates();
		// SecondaryStructFinger s = new SecondaryStructFinger(list.get(0)._2.getCoordinates());

		SecondaryStructFinger s = new SecondaryStructFinger(read("3WST.G"));
		try (Scanner scan = new Scanner(System.in)) {
			String in;
			while (!(in = scan.next()).equals("X")) {
				if (in.equals("g")) {
					int st = scan.nextInt();
					System.out.println(distsToString(s.getRange(st, scan.nextInt()), st));
				}
				else if (in.equals("a"))
					s.printHelices();

			}
		}
		// sc.close();
	}

	public static double[][] dists(Point3d[] pts) {
		double[][] out = new double[pts.length][NUM];
		for (int i = 0; i < NUM; i++) {
			for (int j = 0; j < i; j++) {
				if (pts[i] != null && pts[i - j - 1] != null)
					out[i][j] = pts[i].distance(pts[i - j - 1]);
			}
		}
		for (int i = NUM; i < pts.length; i++) {
			for (int j = 0; j < NUM; j++) {
				if (pts[i] != null && pts[i - j - 1] != null)
					out[i][j] = pts[i].distance(pts[i - j - 1]);
			}
		}
		return out;
	}

	public static String distsToString(double[][] dists) {
		return distsToString(dists, 0);
	}

	public static String distsToString(double[][] dists, int offset) {
		if (dists == null)
			return null;
		StringBuilder a;
		a = new StringBuilder();
		for (int i = 0; i < dists.length; i++)
			a.append((i + offset) + "\t");
		a.append(System.lineSeparator());
		for (int j = 0; j < NUM; j++) {
			for (int i = 0; i < dists.length; i++)
				a.append(String.format("%.3f", dists[i][j]) + "\t");
			a.append(System.lineSeparator());
		}
		return a.toString();
	}

	public static int[] alphaHelices(double[][] dists, int filter) {
		if (filter == 0)
			filter = Integer.MAX_VALUE;
		List<Integer> o = new ArrayList<>();
		boolean on = false;
		for (int i = 0; i < dists.length; i++) {
			if (5 < dists[i][3] && dists[i][3] < 7 && 4 < dists[i][2] && dists[i][2] < 6) {
				if (!on)
					o.add(i);
				on = true;
			}
			else {
				if (on)
					o.add(i);
				on = false;
			}
		}
		int[] out = new int[o.size()];
		for (int i = 0; i < o.size(); i++)
			out[i] = o.get(i);
		return out;
	}

	public static void store(Point3d[] pts, String name) {
		try (PrintWriter pw = new PrintWriter(new FileWriter(name + ".txt"))) {
			pw.println(pts.length);
			for (Point3d p : pts)
				if (p == null)
					pw.println();
				else
					pw.printf("%.3f %.3f %.3f" + System.lineSeparator(), p.x, p.y, p.z);
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static Point3d[] read(String name) {
		Point3d[] pts = null;
		try (Scanner scan = new Scanner(new File(name + ".txt"))) {
			pts = new Point3d[Integer.parseInt(scan.nextLine())];
			for (int i = 0; i < pts.length; i++) {
				String line = scan.nextLine();
				if (line.length() < 2)
					pts[i] = null;
				else {
					String[] spl = line.split(" ");
					pts[i] = new Point3d(Double.parseDouble(spl[0]), Double.parseDouble(spl[1]),
							Double.parseDouble(spl[2]));
				}
			}
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return pts;
	}
}
