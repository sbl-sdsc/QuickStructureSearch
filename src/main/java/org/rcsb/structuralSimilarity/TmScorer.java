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

public class TmScorer {
	private static final String CA_NAME = "CA";
	private static final String GROUP_NAME = "GLU";

	public static Float[] getFatCatTmScore(Point3d[] points1, Point3d[] points2) {
		Float[] scores = new Float[6];
		
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
			return scores;
		}        

		scores[0] = (float) afp.getTMScore();
		scores[1] = (float) afp.getTotalRmsdOpt();
		scores[2] = (float) afp.getProbability();
		scores[3] = (float) afp.getOptLength();
		scores[4] = (float) afp.getCoverage1();
		scores[5] = (float) afp.getCoverage2();
		return scores;
	}

	private static Atom[] getCAAtoms(Point3d[] points) {
		int gaps = 0;
		for (Point3d p: points) {
			if (p == null) {
				gaps++;
			}
		}
		Chain c = new ChainImpl();
		c.setChainID("A");		

		Atom[] atoms = new Atom[points.length-gaps];

		for (int i = 0, j = 0; i < points.length; i++) {
			if (points[i] != null) {
				atoms[j] = new AtomImpl();
				atoms[j].setName(CA_NAME);
				Group g = new AminoAcidImpl();
				g.setPDBName(GROUP_NAME);
				g.addAtom(atoms[j]);
				c.addGroup(g);

				atoms[j].setX(points[i].x);
				atoms[j].setY(points[i].y);
				atoms[j].setZ(points[i].z);
				j++;
			}
		}

		return atoms;
	}
}
