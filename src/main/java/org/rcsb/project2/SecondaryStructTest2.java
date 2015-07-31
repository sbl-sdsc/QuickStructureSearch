package org.rcsb.project2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;
import org.rcsb.structuralSimilarity.ChainSmoother;
import org.rcsb.structuralSimilarity.SavitzkyGolay7PointSmoother;

public class SecondaryStructTest2 {

	private static final String NAME = "2GLP.C";// "2GLL.A";//
	private static final String chainID = NAME.split("\\.")[1];
	private static final boolean smooth = true;
	private static final ChainSmoother CS = new SavitzkyGolay7PointSmoother(1);

	private static final String CA_NAME = "CA";
	private static final String GROUP_NAME = "GLU";

	public static void main(String[] args) {
		Point3d[] pts = SecondaryStructTools.obtain(NAME);
		if (smooth)
			pts = CS.getSmoothedPoints(pts);
		SecondaryStruct s = new SecondaryStruct(pts, smooth);
		File f = new File("data/smooth_" + NAME + ".pdb");
		if (!f.exists()) {
			Chain c = new ChainImpl();
			c.setChainID(NAME.split("\\.")[1]);
			try (PrintWriter pw = new PrintWriter(f)) {
				for (int i = 0; i < pts.length; i++) {
					Atom a = new AtomImpl();
					a = new AtomImpl();
					a.setName(CA_NAME);
					a.setAltLoc(' ');
					Group g = new AminoAcidImpl();
					g.setPDBName(GROUP_NAME);
					g.addAtom(a);
					g.setResidueNumber(chainID, i + 1, null);
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
		// SecondaryStructureSequenceFeature sf = s.getSequenceFeature();
		System.out.println("Start");
		System.out.println(NAME);
		System.out.println(s.getAlphaLength());
		try (Scanner scan = new Scanner(System.in)) {
			String in;
			while (!(in = scan.next()).equals("X")) {
				if (in.equals("g")) {
					int st = scan.nextInt();
					System.out.println(SecondaryStructTools.distsToString(s.getRange(st - 1, scan.nextInt()), st));
				}
				else if (in.equals("a"))
					s.printHelices();
				else if (in.equals("b"))
					s.printStrands();
				else if (in.equals("l"))
					System.out.println(s.length());
				else if (in.equals("c"))
					s.printPoints();
				// else if (in.equals("sf"))
				// for (int i = 0; i < s.length(); i++)
				// System.out.println((i + 1) + ":\t" + sf.toString(i));
				else if (in.equals("test"))
					s.testPrint(s.getAlpha().getFeatures());
				else if (in.equals("test1"))
					System.out.println("=(0,0,0)\t=" + s.normP + "*50\t=" + s.normX + "*50");
				else if (in.equals("test2"))
					SecondaryStruct.printProjection(s.getAlphaNormProjection((byte) 0b00000000));
			}
		}
		// sc.close();
	}
}
