package org.rcsb.project2;

import java.io.File;
import java.io.IOException;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.StructureException;

public class SecondaryStructTestPlotFeat {
	private static String[] chains = new String[] { "1E6R.B", "1E6Z.A" };

	public static void main(String[] args) {
		for (String chain : chains) {
			System.out.println(chain);
			File f = new File("data/" + chain + ".txt");
			SecondaryStruct s = null;
			if (f.exists())
				s = new SecondaryStruct(SecondaryStructTools.read(chain));
			else {
				Point3d[] pts = null;
				try {
					pts = SecondaryStructTools.pull(chain);
				}
				catch (IOException | StructureException e) {
					e.printStackTrace();
				}
				SecondaryStructTools.write(pts, chain);
				s = new SecondaryStruct(pts);
			}
			boolean projection = true;
			boolean alpha = true;
			int i = 0;
			boolean norm = true;
			byte negFlip = (byte) 0b00000000;
			if (projection) {
				if (norm) {
					SecondaryStruct.printProjection(alpha ? s.getAlphaNormProjection(negFlip) : s
							.getBetaNormProjection(negFlip));
				}
				else {
					SecondaryStruct.printProjection(alpha ? s.getAlphaProjection(i) : s.getBetaProjection(i));
				}
			}
			System.out.println();
			s.testPrint(alpha ? s.getHelices() : s.getStrands());
			System.out.println();
			System.out.println("=(0,0,0)\t=" + s.normP + "*50\t=" + s.normX + "*50");
			System.out.println();
			System.out.println();
			System.out.println();
			System.out.println();
		}
	}
}
