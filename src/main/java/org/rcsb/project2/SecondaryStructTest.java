package org.rcsb.project2;

import javax.vecmath.Point3d;

public class SecondaryStructTest {

	public static void main(String[] args) {
		String[] id1 = new String[] { "2VCQ.B", "4FJ4.B", "4FNL.I", "4V6B.Bh", "4V41.P", "4V9N.BV", "3J3Q.1n" };
		String[] id2 = new String[] { "2VCX.B", "4FMW.B", "4FQR.u", "4V6B.CM", "4V45.G", "4VUB.A", "3J3Q.ee" };
		final int N = id1.length;
		int start = 0;
		boolean s = false;
		for (int i = s ? start : 0; i < (s ? (start + 1) : N); i++) {
			Point3d[] p1 = SecondaryStruct.read(id1[i]);
			Point3d[] p2 = SecondaryStruct.read(id2[i]);
			SecondaryStruct s1 = new SecondaryStruct(p1);
			SecondaryStruct s2 = new SecondaryStruct(p2);
			if (s) {
				int i1 = 1;
				int i2 = 1;
				boolean alpha = false;
				if (alpha) {
					System.out.println(id1[i] + "\t\t\t");
					s1.printAlphaProjection(i1);
					System.out.println(id2[i] + "\t\t\t");
					s2.printAlphaProjection(i2);
				}
				else {
					System.out.println(id1[i] + "\t\t\t");
					s1.printBetaProjection(i1);
					System.out.println(id2[i] + "\t\t\t");
					s2.printBetaProjection(i2);
				}
				SecondaryStruct.align(s1, s2, id1[i], id2[i], alpha, i1, i2);
			}
			else {
				SecondaryStruct.align(s1, s2, id1[i], id2[i], true);
				System.out.println();
				SecondaryStruct.align(s1, s2, id1[i], id2[i], false);
				System.out.println();
			}
		}
	}
}
