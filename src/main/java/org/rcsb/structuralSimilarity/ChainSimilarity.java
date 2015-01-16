package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;

public class ChainSimilarity {
	private static final String CA_NAME = "CA";
	private static final String GROUP_NAME = "GLU";

	public static double getFatCatTmScore(Point3d[] points1, Point3d[] points2) {
		Atom[] ca1 = getCAAtoms(points1);
		Atom[] ca2 = getCAAtoms(points2);
		
		FatCatParameters params = new FatCatParameters();
		AFPChain afp = null;
		try {
			StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);
			afp = algorithm.align(ca1,ca2,params);
			double tmScore = AFPChainScorer.getTMScore(afp, ca1, ca2);
			afp.setTMScore(tmScore);
		} catch (StructureException e) {
			e.printStackTrace();
			return 0.0;
		}            
		return afp.getTMScore();
	}

	private static Atom[] getCAAtoms(Point3d[] points) {
		Chain c = new ChainImpl();
		c.setChainID("A");		

		Atom[] atoms = new Atom[points.length];

		for (int i = 0; i < points.length; i++) {
			atoms[i] = new AtomImpl();
			atoms[i].setName(CA_NAME);
			Group g = new AminoAcidImpl();
			g.setPDBName(GROUP_NAME);
			g.addAtom(atoms[i]);
			c.addGroup(g);

			atoms[i].setX(points[i].x);
			atoms[i].setY(points[i].y);
			atoms[i].setZ(points[i].z);
		}

		return atoms;
	}
}
