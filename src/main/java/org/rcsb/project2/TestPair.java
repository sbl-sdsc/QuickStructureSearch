package org.rcsb.project2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import javax.vecmath.Point3d;

import org.rcsb.structuralSimilarity.ChainSmoother;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;

import scala.Tuple2;

/**
 * For testing a pair of proteins
 * 
 * @author Kevin Wu
 *
 */
public class TestPair {

	private static final boolean smooth = true;
	private static final ChainSmoother cs = new SavitzkyGolay7PointSmoother(1);

	public static void main(String[] args) {
		// String[] id1 = new String[] { "2VCQ.B", "4FJ4.B", "4FNL.I", "4V6B.Bh", "4V41.P", "4V9N.BV", "3J3Q.1n" };
		// String[] id2 = new String[] { "2VCX.B", "4FMW.B", "4FQR.u", "4V6B.CM", "4V45.G", "4VUB.A", "3J3Q.ee" };
		String[] id1, id2;
		id1 = id2 = new String[0];
		int N = 0;
		try (BufferedReader br = new BufferedReader(new FileReader("data/pairs.txt"))) {
			N = Integer.parseInt(br.readLine());
			id1 = new String[N];
			id2 = new String[N];
			for (int i = 0; i < N; i++) {
				String[] spl = br.readLine().split(",");
				id1[i] = spl[0];
				id2[i] = spl[1];
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		boolean s = false;
		id1 = new String[] { "4PJX.D" };
		id2 = new String[] { "1IM3.F" };
		N = id1.length;
		s = true;
		int start = 0;
		long t = System.nanoTime();
		for (int i = s ? start : 0; i < (s ? (start + 1) : N); i++) {
			long rt = System.nanoTime();
			Point3d[] p1 = SecondaryStructTools.obtain(id1[i]);
			Point3d[] p2 = SecondaryStructTools.obtain(id2[i]);
			if (smooth) {
				p1 = cs.getSmoothedPoints(p1);
				p2 = cs.getSmoothedPoints(p2);
			}
			SecondaryStruct s1 = new SecondaryStruct(p1, smooth);
			SecondaryStruct s2 = new SecondaryStruct(p2, smooth);
			t += System.nanoTime() - rt;
			if (s) {
				int i1 = 0;
				int i2 = 0;
				boolean alpha = false;
				System.out.println(id1[i] + " vs " + id2[i]);
				if (alpha) {
					// s1.printHelices();
					// s2.printHelices();
					// System.out.println(id1[i] + "\t\t\t");
					// s1.printAlphaProjection(i1);
					// System.out.println(id2[i] + "\t\t\t");
					// s2.printAlphaProjection(i2);
					// SecondaryStructTools.rmsd(s1.getAlphaProjection(i1), s2.getAlphaProjection(i2));
				}
				else {
					// s1.printStrands();
					// s2.printStrands();
					// System.out.println(id1[i] + "\t\t\t");
					// s1.printBetaProjection(i1);
					// System.out.println(id2[i] + "\t\t\t");
					// s2.printBetaProjection(i2);
					// System.out.println();
					// System.out.println();
					Tuple2<Double, int[]> q = SecondaryStructTools.rmsd(s1.getBetaProjection(i1),
							s2.getBetaProjection(i2));
					System.out.println(q._1);
					System.out.println(Arrays.toString(q._2));
					System.out.println("ASDF");
					System.out.println("ASDF");
					System.out.println("ASDF");
				}
				// int al = s1.getAlphaLength() + s2.getAlphaLength();
				int bl = s1.getBetaLength() + s2.getBetaLength();
				// System.out.println("Alen: " + al);
				System.out.println("Blen: " + bl);
				// System.out.println("ScoreA: " + SecondaryStructTools.align(s1.getAlpha(), s2.getAlpha()));
				System.out.println("ScoreB: " + SecondaryStructTools.align(s1.getBeta(), s2.getBeta()));
				// System.out.println("Score: " + SecondaryStructTools.align(s1, s2, alpha));

				// SecondaryStruct.align(s1, s2, id1[i], id2[i], alpha, i1, i2);
			}
			else {
				System.out.println(id1[i] + " vs " + id2[i]);
				// System.out.println("Align: " + SecondaryStructTools.align(s1, s2));
				double a = SecondaryStructTools.align(s1.getAlpha(), s2.getAlpha());
				double b = SecondaryStructTools.align(s1.getBeta(), s2.getBeta());
				// // System.out.println("TM: " + TmScorer.getFatCatTmScore(s1.getPoints(), s2.getPoints())[0]);
				System.out.println("A: " + a);
				System.out.println("B: " + b);
				int al = s1.getAlphaLength() + s2.getAlphaLength();
				int bl = s1.getBetaLength() + s2.getBetaLength();
				System.out.println("Alen: " + al);
				System.out.println("Blen: " + bl);
				double c = (a * al + b * bl) / (al + bl);
				System.out.println("COMB: " + c);
				// System.out.println();
				// System.out.println();
			}
		}
		System.out.println("Time: " + (System.nanoTime() - t) / 1E9);
	}
}
