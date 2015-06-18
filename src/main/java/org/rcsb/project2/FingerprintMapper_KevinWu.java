package org.rcsb.project2;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPChainScorer;

/**
 * FingerprintMapper gives a visual representation of the alignment of two protein chains.
 * 
 * @author Kevin Wu
 */
public final class FingerprintMapper_KevinWu {

	private static final char MATCH_CHAR = '|';
	private static final char MATCH_FRAG_CHAR = '$';
	private static final String SPACER = "\t";
	private static final String CA_NAME = "CA";
	private static final String GROUP_NAME = "GLU";

	private FingerprintMapper_KevinWu() {
	}

	/**
	 * Gets the Optimal Alignment of 2 protein chains, including gaps
	 * 
	 * @param points1
	 *            List of Point3d representing the first protein chain.
	 * @param points2
	 *            List of Point3d representing the second protein chain.
	 * @return 3 dimensional array for the optimal alignment. Includes gaps.
	 */
	private static int[][][] getOptimalAlignment(Point3d[] points1, Point3d[] points2) {
		List<Integer> gaps1 = new ArrayList<>();
		List<Integer> gaps2 = new ArrayList<>();
		Atom[] ca1 = getCAAtoms(points1, gaps1);
		Atom[] ca2 = getCAAtoms(points2, gaps2);

		FatCatParameters params = new FatCatParameters();
		AFPChain afp = null;
		try {
			StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);
			afp = algorithm.align(ca1, ca2, params);
			double tmScore = AFPChainScorer.getTMScore(afp, ca1, ca2);
			afp.setTMScore(tmScore);
		}
		catch (StructureException e) {
			e.printStackTrace();
		}
		if (afp == null) {
			new NullPointerException("AFP is null").printStackTrace();
			return new int[0][0][0];
		}
		int[][][] withGaps = afp.getOptAln();
		addGaps(withGaps, gaps1, gaps2, 0);
		return withGaps;
	}

	/**
	 * Adds gaps into the Optimal Alignment given by FATCAT.<br>
	 * FATCAT can give multiple alignments, this inserts the gaps into the one given by <code>ind</code>.
	 * 
	 * @param optAln
	 *            3 dimensional Array of the optimal alignment given by
	 *            {@link org.biojava.nbio.structure.align.model.AFPChain#getOptAln()}
	 * @param gaps1
	 *            List of Integer indexes of gaps in protein chain 1
	 * @param gaps2
	 *            List of Integer indexes of gaps in protein chain 2
	 * @param index
	 *            index in optAln to insert gaps
	 */
	private static void addGaps(int[][][] optAln, List<Integer> gaps1, List<Integer> gaps2, int index) {
		final int N = Math.max(gaps1.size(), gaps2.size()) + optAln[0][0].length;
		gaps1.add(N);
		gaps2.add(N);
		int[][] wg = new int[2][N];
		for (int i = 0, j = 0; i < gaps1.size(); i++) {
			for (; j < gaps1.get(i); j++)
				wg[0][j] = optAln[index][0][j - i];
			if (j != N)
				wg[0][j++] = -1;
		}
		for (int i = 0, j = 0; i < gaps2.size(); i++) {
			for (; j < gaps2.get(i); j++)
				wg[1][j] = optAln[index][1][j - i];
			if (j != N)
				wg[1][j++] = -1;
		}
		optAln[index] = wg;
	}

	/**
	 * Gets Array of Atom for FATCAT algorithm. <br>
	 * Also records the indexes of gaps in a gaps List.
	 * 
	 * @param points
	 *            Array of Point3d for the points in the protein chain.
	 * @param gaps
	 *            List to add gap indexes to.
	 * @return Array of Atom for FATCAT algorithm.
	 */
	private static Atom[] getCAAtoms(Point3d[] points, List<Integer> gaps) {
		if (gaps.size() != 0) {
			new IllegalArgumentException(String.format(
					"Error: gaps (size %d) already contains some elements. It will be cleared.", gaps.size()))
					.printStackTrace();
			gaps.clear();
		}

		for (int i = 0; i < points.length; i++)
			if (points[i] == null)
				gaps.add(i);

		Chain c = new ChainImpl();
		c.setChainID("A");

		Atom[] atoms = new Atom[points.length - gaps.size()];

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

	/**
	 * Gives a String that shows the alignment of two protein chains with given fragments. <br>
	 * <br>
	 * Row 1 has fragment IDs for chain 1<br>
	 * Row 2 has positions of points in chain 1<br>
	 * Row 3 has matching<br>
	 * {@value #MATCH_CHAR} means the points at the position match (based on fatcat) <br>
	 * {@value #MATCH_FRAG_CHAR} means the fragments at the position match<br>
	 * Row 4 has positions of points in chain 2<br>
	 * Row 5 has fragment IDs for chain 2<br>
	 *
	 * @param points1
	 *            Array of Point3d representing the points of the first protein chain
	 * @param points2
	 *            Array of Point3d representing the points of the second protein chain
	 * @param seqFeat1
	 *            SequenceFeatureInterface for the first protein chain
	 * @param seqFeat2
	 *            SequenceFeatureInterface for the second protein chain
	 * @return String showing the alignment
	 */
	public static <T> String align(Point3d[] points1, Point3d[] points2, SequenceFeatureInterface<T> seqFeat1,
			SequenceFeatureInterface<T> seqFeat2) {
		int[][][] optAln = getOptimalAlignment(points1, points2);
		int[] topMatchIndexes = optAln[0][0];
		int[] botMatchIndexes = optAln[0][1];

		StringBuilder topFrag = new StringBuilder();
		StringBuilder topIndexes = new StringBuilder();
		StringBuilder matches = new StringBuilder();
		StringBuilder botIndexes = new StringBuilder();
		StringBuilder botFeats = new StringBuilder();
		for (int matchInd = 0; matchInd < topMatchIndexes.length; matchInd++) {
			if (topMatchIndexes[matchInd] != 0 || botMatchIndexes[matchInd] != 0 || matchInd == 0) {
				matches.append(MATCH_CHAR);
				topIndexes.append(topMatchIndexes[matchInd]);
				botIndexes.append(botMatchIndexes[matchInd]);
			}
			else {
				if (topMatchIndexes[matchInd] == -1)
					topIndexes.append("nul");
				else
					topIndexes.append(String.format("(%d)", matchInd + topMatchIndexes[0]));
				if (botMatchIndexes[matchInd] == -1)
					botIndexes.append("nul");
				else
					botIndexes.append(String.format("(%d)", matchInd + botMatchIndexes[0]));
			}
			if (seqFeat1 != null && seqFeat2 != null)
				if (seqFeat1.identity(seqFeat2, matchInd, matchInd))
					matches.append(MATCH_FRAG_CHAR);
				else
					matches.append(String.format("%.2f", seqFeat1.similarity(seqFeat2, matchInd, matchInd)));
			topFrag.append(seqFeat1 == null ? "" : seqFeat1.get(matchInd));
			botFeats.append(seqFeat2 == null ? "" : seqFeat2.get(matchInd));
			matches.append(SPACER);
			topIndexes.append(SPACER);
			botIndexes.append(SPACER);
			topFrag.append(seqFeat1 == null ? "" : SPACER);
			botFeats.append(seqFeat1 == null ? "" : SPACER);
		}
		topFrag.append(seqFeat1 == null ? "" : System.lineSeparator());
		topFrag.append(topIndexes);
		topFrag.append(System.lineSeparator());
		topFrag.append(matches);
		topFrag.append(System.lineSeparator());
		topFrag.append(botIndexes);
		topFrag.append(seqFeat2 == null ? "" : System.lineSeparator());
		topFrag.append(seqFeat2 == null ? "" : botFeats);
		return topFrag.toString();
	}

	/**
	 * Gives a String that shows the alignment of two protein chains (no fragments). <br>
	 * <br>
	 * Row 1 has positions of points in chain 1<br>
	 * Row 2 has matching<br>
	 * {@value #MATCH_CHAR} means the points at the position match (based on fatcat) <br>
	 * Row 3 has positions of points in chain 2<br>
	 * 
	 * @param points1
	 *            Array of Point3d representing the points of the first protein chain
	 * @param points2
	 *            Array of Point3d representing the points of the second protein chain
	 * @return String showing the alignment
	 */
	public static String align(Point3d[] points1, Point3d[] points2) {
		return align(points1, points2, null, null);
	}
}
