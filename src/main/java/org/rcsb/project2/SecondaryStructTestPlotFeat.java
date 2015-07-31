package org.rcsb.project2;

import javax.vecmath.Point3d;

import org.rcsb.structuralSimilarity.ChainSmoother;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;

public class SecondaryStructTestPlotFeat {
	private static String[] chains = new String[] { "3VKE.A", "1PUE.E" };
	private static final boolean smooth = true;
	private static final boolean write = true;
	private static final ChainSmoother CS = new SavitzkyGolay7PointSmoother(1);

	public static void main(String[] args) {
		for (String chain : chains) {
			Point3d[] pts = SecondaryStructTools.obtain(chain);
			System.out.println(chain + "\t\t\t\t");
			if (smooth)
				pts = CS.getSmoothedPoints(pts);
			SecondaryStruct s = new SecondaryStruct(pts, smooth);
			boolean projection = false;
			boolean alpha = false;
			boolean plotPts = false;
			// int i = 0;
			boolean norm = true;
			byte negFlip = (byte) 0b00000000;
			if (projection) {
				if (norm) {
					SecondaryStruct.printProjection(alpha ? s.getAlphaNormProjection(negFlip) : s
							.getBetaNormProjection(negFlip));
				}
				else {
					// SecondaryStruct.printProjection(alpha ? s.getAlphaProjection(i) : s.getBetaProjection(i));
				}
			}
			if (plotPts) {
				for (int j = 0; j < pts.length; j++) {
					System.out.printf("=%s\t=A%d-A%d\t=Segment[B%d,B%d]" + System.lineSeparator(), pts[j], j + 2,
							pts.length + 2, j, j + 1);
				}
				System.out.printf("=Sum[A2:A%d]/%d", pts.length + 1, pts.length);
			}
			if (write) {
			}
			System.out.println();
			s.testPrint(alpha ? s.getAlpha().getFeatures() : s.getBeta().getFeatures());
			System.out.println();
			System.out.println("=(0,0,0)\t=" + s.normP + "*50\t=" + s.normX + "*50");
			System.out.println();
			System.out.println();
			System.out.println();
			System.out.println();
		}
	}
}
