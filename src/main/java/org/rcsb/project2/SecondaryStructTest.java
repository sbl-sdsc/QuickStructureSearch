package org.rcsb.project2;

import javax.vecmath.Point3d;

public class SecondaryStructTest {

	public static void main(String[] args) {
		String[] id1 = new String[] { "2VCQ.B", "4FJ4.B", "4FNL.I", "4V6B.Bh", "4V41.P", "4V9N.BV", "3J3Q.1n" };
		String[] id2 = new String[] { "2VCX.B", "4FMW.B", "4FQR.u", "4V6B.CM", "4V45.G", "4VUB.A", "3J3Q.ee" };
		// final int N = id1.length;
		int start = 2;
		for (int i = start; i < start + 1; i++) {
			Point3d[] p1 = SecondaryStructFinger.read(id1[i]);
			Point3d[] p2 = SecondaryStructFinger.read(id2[i]);
			SecondaryStructFinger s1 = new SecondaryStructFinger(p1);
			SecondaryStructFinger s2 = new SecondaryStructFinger(p2);
			// int nA = Math.min(s1.getAlphaLength(), s2.getAlphaLength());
			if (s1.getAlphaLength() != 0) {
				System.out.println("Alpha Helices");
				int d1, d2;
				d1 = d2 = 0;
				for (int k = 0; k < 1 /* nA */; k++) {
					System.out.println((k + d1) + " vs " + (k + d2));
					int ind1, ind2;
					ind1 = ind2 = k;
					ind1 += d1;
					ind2 += d2;
					System.out.println(id1[i] + "\t\t\t");
					s1.printAlphaProjection(ind1);
					System.out.println(id2[i] + "\t\t\t");
					s2.printAlphaProjection(ind2);
					System.out.println("RMSD: "
							+ SecondaryStructFinger.rmsd(s1.getAlphaProjection(ind1), s2.getAlphaProjection(ind2)));
					// System.out.println(SecondaryStructFinger.rmsd(s1.getAlphaProjection(ind1),
					// s2.getAlphaProjection(ind2)));
					// SecondaryStructProjection a1 = s1.getAlphaProjection(ind1);
					// SecondaryStructProjection a2 = s2.getAlphaProjection(ind2);
					// a1.getCloseTo(a2.getStart(0), a2.getEnd(0));
					System.out.println();
				}
			}
			else {
				System.out.println("no alpha helices");
			}
			if (s1.getBetaLength() != 0 && s2.getBetaLength() != 0) {
				System.out.println("Beta Strands");
				System.out.println(s1.getBetaLength());
				System.out.println(s2.getBetaLength());
				int d1, d2;
				d1 = d2 = 0;
				for (int k = 0; k < 1/* nA */; k++) {
					System.out.println((k + d1) + " vs " + (k + d2));
					int ind1, ind2;
					ind1 = ind2 = k;
					ind1 += d1;
					ind2 += d2;
					System.out.println(id1[i] + "\t\t\t");
					s1.printBetaProjection(ind1);
					System.out.println(id2[i] + "\t\t\t");
					s2.printBetaProjection(ind2);
					System.out.println("RMSD: "
							+ SecondaryStructFinger.rmsd(s1.getBetaProjection(ind1), s2.getBetaProjection(ind2)));
					// System.out.println(SecondaryStructFinger.rmsd(s1.getBetaProjection(ind1),
					// s2.getBetaProjection(ind2)));
					// SecondaryStructProjection a1 = s1.getBetaProjection(ind1);
					// SecondaryStructProjection a2 = s2.getBetaProjection(ind2);
					// a1.getCloseTo(a2.getStart(0), a2.getEnd(0));
					System.out.println();
				}

			}
			else {
				System.out.println("no beta strands");
			}
			System.out.println();
		}
	}

}
