package org.rcsb.project2;

import java.util.Arrays;

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
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.project3.SequenceFeatureInterface;

/**
 * FingerprintMapper gives a visual representation of the alignment of two protein chains.
 * 
 * @author Kevin Wu
 */
public final class FingerprintMapper_KevinWu {
	private static final SequenceFeatureInterface<String> BLANK_SEQUENCE_FEATURE = new SequenceFeatureInterface<String>() {

		@Override
		public double similarity(SequenceFeatureInterface<String> sequence2, int i, int j) {
			return 0;
		}

		@Override
		public boolean identity(SequenceFeatureInterface<String> sequence2, int i, int j) {
			return true;
		}

		@Override
		public String[] getSequence() {
			return null;
		}

		@Override
		public String get(int index) {
			return "";
		}

		@Override
		public int length() {
			return 0;
		}

		@Override
		public String toString(int index) {
			return "";
		}

		@Override
		public double todouble(int index) {
			return 0;
		}

	};
	private static final char MATCH_CHAR = '|';
	private static final char SIMILAR = ':';
	private static final char NOT_SIMILAR = '.';
	private static final char MATCH_FRAG_CHAR = '$';
	private static final String SPACER = "\t";
	private static final String CA_NAME = "CA";
	private static final String GROUP_NAME = "GLU";
	private static final double MATCH_THRESHOLD = 0.5;

	private FingerprintMapper_KevinWu() {
	}

	/**
	 * Gets the Optimal Alignment of 2 protein chains, including gaps
	 * 
	 * @param points1
	 *            List of Point3d for protein chain 1.
	 * @param points2
	 *            List of Point3d for protein chain 2.
	 * @return 3 dimensional array for the optimal alignment. Includes gaps.
	 */
	private static int[][][] getOptimalAlignment(Point3d[] points1, Point3d[] points2) {
		Atom[] ca1 = getCAAtoms(points1);
		Atom[] ca2 = getCAAtoms(points2);
		FatCatParameters params = new FatCatParameters();
		AFPChain afp = null;
		try {
			StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);
			afp = algorithm.align(ca1, ca2, params);
			// double tmScore = AFPChainScorer.getTMScore(afp, ca1, ca2);
			// afp.setTMScore(tmScore);
		}
		catch (StructureException e) {
			e.printStackTrace();
		}
		if (afp == null) {
			new NullPointerException("AFP is null").printStackTrace();
			return new int[0][0][0];
		}
		int[][][] opt = afp.getOptAln();
		return opt;
	}

	/**
	 * Gets Array of Atom for FATCAT algorithm. <br>
	 * 
	 * @param points
	 *            Array of Point3d for the protein chain.
	 * @return Array of Atom for FATCAT algorithm.
	 */
	private static Atom[] getCAAtoms(Point3d[] points) {
		int gaps = 0;
		for (int i = 0; i < points.length; i++)
			if (points[i] == null)
				gaps++;

		Chain c = new ChainImpl();
		c.setChainID("A");

		Atom[] atoms = new Atom[points.length - gaps];

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
	 * Returns String showing the alignment of 2 protein chains using FATCAT and given {@link SequenceFeatureInterface}s <br>
	 * <br>
	 * Row 1 has fragment IDs for chain 1<br>
	 * Row 2 has positions of points in chain 1<br>
	 * Row 3 has matching<br>
	 * {@value #MATCH_CHAR} means the points at the position match (based on fatcat) <br>
	 * {@value #MATCH_FRAG_CHAR} means the fragments at the position match<br>
	 * Row 4 has positions of points in chain 2<br>
	 * Row 5 has fragment IDs for chain 2<br>
	 *
	 * @param p1
	 *            Array of Point3d for protein chain 1
	 * @param p2
	 *            Array of Point3d for protein chain 2
	 * @param f1
	 *            {@link SequenceFeatureInterface} for protein chain 1.
	 * @param f2
	 *            {@link SequenceFeatureInterface} for protein chain 2.
	 * @param letter
	 *            true for 1 letter, false for verbose
	 * @return String showing the alignment
	 */
	public static <T> String align(Point3d[] p1, Point3d[] p2, SequenceFeatureInterface<T> f1,
			SequenceFeatureInterface<T> f2, boolean letter) {
		int[][][] optAln = getOptimalAlignment(p1, p2);
		int[] a1, a2;
		if (optAln != null) {
			a1 = optAln[0][0];
			a2 = optAln[0][1];
		}
		else {
			a1 = new int[0];
			a2 = new int[0];
		}
		System.out.println((Arrays.toString(a1) + System.lineSeparator() + Arrays.toString(a2)).replaceAll(" ", "\t"));
		final String[] n1 = new String[p1.length];
		for (int i = 0; i < n1.length; i++)
			n1[i] = Integer.toString(i);
		return align(a1, a2, n1, n1, f1, f2, letter);
	}

	/**
	 * Gives a String that shows the alignment of two protein chains (no fragments). <br>
	 * <br>
	 * Row 1 has positions of points in chain 1<br>
	 * Row 2 has matching<br>
	 * {@value #MATCH_CHAR} means the points at the position match (based on fatcat) <br>
	 * Row 3 has positions of points in chain 2<br>
	 * 
	 * @param p1
	 *            Array of Point3d for protein chain 1.
	 * @param p2
	 *            Array of Point3d for protein chain 2.
	 * @param letter
	 *            true for 1 letter, false for verbose
	 * @return String showing the alignment
	 */
	public static String align(Point3d[] p1, Point3d[] p2, boolean letter) {
		return align(p1, p2, BLANK_SEQUENCE_FEATURE, BLANK_SEQUENCE_FEATURE, letter);
	}

	/**
	 * Returns String showing the alignment of 2 protein chains using FATCAT and given {@link SequenceFeatureInterface}s
	 * 
	 * @param c1
	 *            {@link SimplePolymerChain} for protein chain 1.
	 * @param c2
	 *            {@link SimplePolymerChain} for protein chain 2.
	 * @param f1
	 *            {@link SequenceFeatureInterface} for protein chain 1.
	 * @param f2
	 *            {@link SequenceFeatureInterface} for protein chain 2.
	 * @param letter
	 *            true for 1 letter, false for verbose
	 * @return String showing the alignment.
	 */
	public static <T> String align(SimplePolymerChain c1, SimplePolymerChain c2, SequenceFeatureInterface<T> f1,
			SequenceFeatureInterface<T> f2, boolean letter) {
		int[][][] optAln = getOptimalAlignment(c1.getCoordinates(), c2.getCoordinates());
		int[] t, b;
		t = b = new int[0];
		if (optAln != null) {
			t = optAln[0][0];
			b = optAln[0][1];
		}
		System.out.println((Arrays.toString(t) + System.lineSeparator() + Arrays.toString(b)).replaceAll(" ", "\t"));
		return align(t, b, c1.getSequence().split(""), c2.getSequence().split(""), f1, f2, letter);
	}

	/**
	 * Returns String showing the alignment of 2 protein chains using FATCAT and given {@link SequenceFeatureInterface}s
	 * 
	 * @param a1
	 *            Position of matching points in protein chain 1.
	 * @param a2
	 *            Position of matching points in protein chain 2.
	 * @param n1
	 *            String for points in protein chain 1.
	 * @param n2
	 *            String for points in protein chain 2.
	 * @param f1
	 *            {@link SequenceFeatureInterface} for protein chain 1.
	 * @param f2
	 *            {@link SequenceFeatureInterface} for protein chain 2.
	 * @param letter
	 *            true for 1 letter, false for verbose
	 * @return String showing the alignment of 2 protein chains.
	 */
	public static <T> String align(int[] a1, int[] a2, String[] n1, String[] n2, SequenceFeatureInterface<T> f1,
			SequenceFeatureInterface<T> f2, boolean letter) {
		String SPACER = letter ? "" : FingerprintMapper_KevinWu.SPACER;
		if (a1.length != a2.length)
			new IllegalArgumentException("alignment lengths not equal (" + a1.length + " ," + a2.length + ")")
					.printStackTrace();
		StringBuilder tf = new StringBuilder(), ti = new StringBuilder(), m = new StringBuilder(), bf = new StringBuilder(), bi = new StringBuilder();
		if (a1.length == 0) {
			noMatch(n1, n2, f1, f2, tf, ti, m, bi, bf, 0, n1.length, 0, n2.length, letter);
		}
		else {

			int tp, bp;
			tp = bp = -1;
			for (int i = 0; i < a1.length; i++) {
				if (i != 0 && (a1[i] == 0 || a2[i] == 0)) {
					noMatch(n1, n2, f1, f2, tf, ti, m, bi, bf, tp + 1, n1.length, bp + 1, n2.length, letter);
					break;
				}
				if (a1[i] - tp > 1 || a2[i] - bp > 1)
					noMatch(n1, n2, f1, f2, tf, ti, m, bi, bf, tp + 1, a1[i], bp + 1, a2[i], letter);
				tf.append(letter ? f1.toString(a1[i]).charAt(0) : f1.toString(a1[i]));
				ti.append(letter ? n1[a1[i]].charAt(0) : n1[a1[i]]);
				if (!letter) {
					m.append(MATCH_CHAR);
					m.append(f1.identity(f2, a1[i], a2[i]) ? MATCH_FRAG_CHAR : String.format("%.2f",
							f1.similarity(f2, a1[i], a2[i])));
				}
				else {
					m.append(f1.identity(f2, a1[i], a2[i]) ? MATCH_FRAG_CHAR
							: f1.similarity(f2, a1[i], a2[i]) > MATCH_THRESHOLD ? SIMILAR : NOT_SIMILAR);
				}
				bi.append(letter ? n2[a2[i]].charAt(0) : n2[a2[i]]);
				bf.append(letter ? f2.toString(a2[i]).charAt(0) : f2.toString(a2[i]));

				tf.append(SPACER);
				ti.append(SPACER);
				m.append(SPACER);
				bi.append(SPACER);
				bf.append(SPACER);
				tp = a1[i];
				bp = a2[i];
			}
		}
		tf.append(System.lineSeparator());
		tf.append(ti);
		tf.append(System.lineSeparator());
		tf.append(m);
		tf.append(System.lineSeparator());
		tf.append(bi);
		tf.append(System.lineSeparator());
		tf.append(bf);
		return tf.toString();
	}

	/**
	 * For gaps
	 * 
	 * @param n1
	 *            String array for points on protein chain 1.
	 * @param n2
	 *            String array for points on protein chain 2.
	 * @param f1
	 *            {@link SequenceFeatureInterface} for protein chain 1.
	 * @param f2
	 *            {@link SequenceFeatureInterface} for protein chain 2.
	 * @param tf
	 *            {@link StringBuilder} for protein 1 fragments.
	 * @param ti
	 *            {@link StringBuilder} for protein 1 chain.
	 * @param m
	 *            {@link StringBuilder} for matches.
	 * @param bi
	 *            {@link StringBuilder} for protein 2 chain.
	 * @param bf
	 *            {@link StringBuilder} for protein 2 fragments.
	 * @param tS
	 *            Protein 1 gap start position (inclusive).
	 * @param tE
	 *            Protein 1 gap end position (exclusive).
	 * @param bS
	 *            Protein 2 gap start position (inclusive).
	 * @param bE
	 *            Protein 2 gap end position (exclusive).
	 * @param letter
	 *            true for 1 letter, false for verbose
	 */
	private static <T> void noMatch(String[] n1, String[] n2, SequenceFeatureInterface<T> f1,
			SequenceFeatureInterface<T> f2, StringBuilder tf, StringBuilder ti, StringBuilder m, StringBuilder bi,
			StringBuilder bf, int tS, int tE, int bS, int bE, boolean letter) {
		String SPACER = letter ? " " : FingerprintMapper_KevinWu.SPACER;
		int tv = tE - tS;
		int bv = bE - bS;
		for (int i = 0; i < tv; i++) {
			ti.append(letter ? n1[i + tS].charAt(0) : n1[i + tS]);
			tf.append(letter ? f1.toString(i).charAt(0) : f1.toString(i));
			if (!letter) {
				ti.append(SPACER);
				tf.append(SPACER);
			}
		}
		for (int i = 0; i < bv; i++) {
			bi.append(letter ? n2[i + bS].charAt(0) : n2[i + bS]);
			bf.append(letter ? f2.toString(i).charAt(0) : f2.toString(i));
			if (!letter) {
				bi.append(SPACER);
				bf.append(SPACER);
			}
		}
		int mv = Math.max(tv, bv);
		for (int i = 0; i < mv; i++)
			// Match spacing
			m.append(SPACER);
		for (int i = 0; i < mv - tv; i++) { // top spacing
			tf.append(SPACER);
			ti.append(SPACER);
		}
		for (int i = 0; i < mv - bv; i++) { // bot spacing
			bf.append(SPACER);
			bi.append(SPACER);
		}
	}
}
