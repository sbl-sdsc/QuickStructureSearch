package org.rcsb.project2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;
import org.rcsb.structuralSimilarity.ChainSmoother;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;

/**
 * For plotting the features of any number of given proteins
 * 
 * @author Kevin Wu
 *
 */
public class TestPlotFeatures {
	private static String[] chains = new String[] { "4V9N.BV", "4VUB.A" };
	private static final boolean smooth = false;
	private static final ChainSmoother CS = new SavitzkyGolay7PointSmoother(2);
	private static final String GROUP_NAME = "GLU";
	private static final String CA_NAME = "CA";

	public static void main(String[] args) {
		for (String chain : chains) {
			Point3d[] pts = SecondaryStructTools.obtain(chain);
			System.out.println(chain + "\t\t\t\t");
			if (smooth)
				pts = CS.getSmoothedPoints(pts);
			SecondaryStruct s = new SecondaryStruct(pts, smooth);
			boolean projection = true;
			boolean alpha = false;
			boolean plotPts = false;
			boolean norm = true;
			byte negFlip = (byte) 0b000000000;
			boolean printFeat = false;
			// int i = 0;
			if (projection) {
				if (norm) {
					SecondaryStruct.printProjection(alpha ? s.getAlpha().getNormProjection(negFlip) : s.getBeta()
							.getNormProjection(negFlip));
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
			System.out.println();
			if (printFeat)
				testPrintFeature(s, chain, false);
			System.out.println();
			// System.out.println("=(0,0,0)\t=" + s.normP + "*50\t=" + s.normX + "*50");
			// System.out.println();
			// System.out.println();
			// System.out.println();
			// System.out.println();
		}
	}

	public static void testPrint(Point3d[] pts, String name) {
		File f = new File("data/smooth_" + name + ".pdb");
		Chain c = new ChainImpl();
		c.setChainID(name.split("\\.")[1]);
		try (PrintWriter pw = new PrintWriter(f)) {
			for (int i = 0; i < pts.length; i++) {
				Atom a = new AtomImpl();
				a = new AtomImpl();
				a.setName(CA_NAME);
				a.setAltLoc(' ');
				Group g = new AminoAcidImpl();
				g.setPDBName(GROUP_NAME);
				g.addAtom(a);
				g.setResidueNumber(name.split("\\.")[1], i + 1, null);
				c.addGroup(g);
				a.setX(pts[i].x);
				a.setY(pts[i].y);
				a.setZ(pts[i].z);
				pw.print(a.toPDB());
			}
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private static void testPrintFeature(SecondaryStruct s, String name, boolean alpha) {
		Point3d[] pts = s.getPoints();
		File f = new File("data/" + name + (alpha ? "_ALPHA" : "_BETA") + ".pdb");
		SecondaryStructFeature ssf = alpha ? s.getAlpha() : s.getBeta();
		int[] st, ed;
		st = ssf.getFeatures()._1;
		ed = ssf.getFeatures()._2;
		try (PrintWriter pw = new PrintWriter(new FileWriter(f))) {
			for (int i = 0; i < ssf.length(); i++) {
				Chain c = new ChainImpl();
				String chainID = String.valueOf((char) ('A' + i));
				c.setChainID(chainID);
				Atom b = new AtomImpl();
				b.setName(CA_NAME);
				b.setAltLoc(' ');
				Group gb = new AminoAcidImpl();
				gb.setPDBName(GROUP_NAME);
				gb.addAtom(b);
				gb.setResidueNumber(chainID, i * 2, null);
				c.addGroup(gb);
				b.setX(pts[st[i]].x);
				b.setY(pts[st[i]].y);
				b.setZ(pts[st[i]].z);
				pw.print(b.toPDB());
				Atom e = new AtomImpl();
				e.setName(CA_NAME);
				e.setAltLoc(' ');
				Group ge = new AminoAcidImpl();
				ge.setPDBName(GROUP_NAME);
				ge.addAtom(e);
				ge.setResidueNumber(chainID, i * 2 + 1, null);
				c.addGroup(ge);
				e.setX(pts[ed[i]].x);
				e.setY(pts[ed[i]].y);
				e.setZ(pts[ed[i]].z);
				pw.print(e.toPDB());
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

	}
}
