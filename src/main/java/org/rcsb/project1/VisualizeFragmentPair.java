package org.rcsb.project1;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;

public class VisualizeFragmentPair {
	private static final String CA_NAME = "CA";
	private static final String GROUP_NAME = "GLU";
	private static final Character iCode = ' ';

	public static void writeFragmentPair(Point3d[] x, Point3d[] y, String fileName) {
		PrintWriter writer = null;
		try {
			writer = new PrintWriter(fileName);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Atom[] ax = getCAAtoms(x, "A");
		Atom[] ay = getCAAtoms(x, "B");
		
		for (Atom atom: ax) {
			writer.println(atom.toPDB());
		}
		for (Atom atom: ay) {
			writer.println(atom.toPDB());
		}
		writer.close();
	}

	private static Atom[] getCAAtoms(Point3d[] points, String chainId) {
		int gaps = 0;
		for (Point3d p: points) {
			if (p == null) {
				gaps++;
			}
		}
		Chain c = new ChainImpl();
		c.setChainID(chainId);		

		Atom[] atoms = new Atom[points.length-gaps];

		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				atoms[i] = new AtomImpl();
				atoms[i].setName(CA_NAME);
				Group g = new AminoAcidImpl();
				g.setPDBName(GROUP_NAME);
				g.setResidueNumber(chainId, i, iCode);
				g.addAtom(atoms[i]);
				c.addGroup(g);

				atoms[i].setX(points[i].x);
				atoms[i].setY(points[i].y);
				atoms[i].setZ(points[i].z);
			}
		}

		return atoms;
	}
}
