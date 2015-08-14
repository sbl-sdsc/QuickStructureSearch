package org.rcsb.project2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.symmetry.geometry.MomentsOfInertia;

import scala.Tuple2;

public class SecondaryStruct {
	private static final int BETA_FILTER = 3;
	private static final int ALPHA_FILTER = 4;
	public static final double C = 3;
	public static final double C_SQ = C * C;

	static final int NUM = 4;

	private Point3d[] pts;
	private boolean smooth;
	private double[][] dists;
	Vector3d normP, normX;
	Point3d normC;
	SecondaryStructFeature alphaHelices, betaStrands;

	protected SecondaryStruct() {
	}

	public SecondaryStruct(Point3d[] pts, boolean smooth) {
		this.smooth = smooth;
		List<Point3d> p = new ArrayList<>();
		for (Point3d pt : pts)
			if (pt != null)
				p.add(pt);
		this.pts = new Point3d[p.size()];
		p.toArray(this.pts);
		dists = SecondaryStructTools.dists(this.pts);
		initNormProjections();
		// System.out.println("alpha");
		alphaHelices = new AlphaBetaStruct(SecondaryStructTools.alphaHelices(dists, this.pts, ALPHA_FILTER, smooth),
				normP, normX, normC, pts);
		// System.out.println("beta");
		betaStrands = new AlphaBetaStruct(SecondaryStructTools.betaStrands(dists, this.pts, BETA_FILTER, smooth),
				normP, normX, normC, pts);
	}

	/**
	 * The number of points in this protein chain.
	 * 
	 * @return The number of points in this protein chain.
	 */
	public int length() {
		return dists.length;
	}

	/**
	 * @return length of alpha helices
	 */
	public int getAlphaLength() {
		return alphaHelices.length();
	}

	/**
	 * @return length of beta strands
	 */
	public int getBetaLength() {
		return betaStrands.length();
	}

	/**
	 * @return number of alpha helix points
	 */
	public int getNumAlphaPoints() {
		return alphaHelices.getNumPoints();
	}

	/**
	 * @return number of beta strand points
	 */
	public int getNumBetaPoints() {
		return betaStrands.getNumPoints();
	}

	/**
	 * Returns the ith set of distances
	 * 
	 * @param i
	 *            index to get the distances from
	 * @return The ith set of distances
	 */
	public double[] get(int i) {
		return dists[i];
	}

	/**
	 * Returns an array with the distances between 2 indeces, or <br>
	 * null if the indeces are bad.
	 * 
	 * @param i
	 *            Starting index (inclusive).
	 * @param j
	 *            Ending index (exclusive).
	 * @return An array with the distances between 2 indeces
	 */
	public double[][] getRange(int i, int j) {
		if (i < 0 || i > dists.length || j < 0 || j > dists.length || i > j)
			return null;
		return Arrays.copyOfRange(dists, i, j);
	}

	/**
	 * @return Point3d array (backbone)
	 */
	public Point3d[] getPoints() {
		return pts;
	}

	/**
	 * @return {@link SecondaryStructFeature} for the alpha helices
	 */
	public SecondaryStructFeature getAlpha() {
		return alphaHelices;
	}

	/**
	 * @return {@link SecondaryStructFeature} for the beta strands
	 */
	public SecondaryStructFeature getBeta() {
		return betaStrands;
	}

	/**
	 * Return a {@link SecondaryStructProjection} for the projection of all the alpha helix vectors onto the ith alpha
	 * helix plane
	 * 
	 * @param i
	 *            index of alpha helix to project onto.
	 * @return A {@link SecondaryStructProjection} for the projection of all the alpha helix vectors onto the ith alpha
	 *         helix plane
	 */
	public SecondaryStructProjection getAlphaProjection(int i) {
		return alphaHelices.getProjection(i);
	}

	/**
	 * Return a {@link SecondaryStructProjection} for the projection of all the beta strand vectors onto the ith beta
	 * strand plane
	 * 
	 * @param i
	 *            index of beta strand to project onto.
	 * @return A {@link SecondaryStructProjection} for the projection of all the beta strand vectors onto the ith beta
	 *         strand plane
	 */
	public SecondaryStructProjection getBetaProjection(int i) {
		return betaStrands.getProjection(i);
	}

	/**
	 * initialize the principal axes vectors
	 */
	private void initNormProjections() {
		MomentsOfInertia m = new MomentsOfInertia();
		for (Point3d pt : pts) {
			m.addPoint(pt, 1.0);
		}
		Point3d c = m.centerOfMass();
		Vector3d[] v = m.getPrincipalAxes();
		double[] vals = m.getPrincipalMomentsOfInertia();
		// System.out.println(Arrays.toString(vals));
		int hi, lo;
		hi = lo = 0;
		for (int i = 0; i < 3; i++) {
			if (vals[i] > hi)
				hi = i;
			if (vals[i] < lo)
				lo = i;
		}
		// System.out.println("hi: " + hi);
		// System.out.println("lo: " + lo);
		normC = new Point3d(c);
		normP = new Vector3d(v[lo]);
		normX = new Vector3d(v[hi]);
	}

	/**
	 * returns the norm projections for the alpha helices using the bits in b to determine which axes to flip
	 * 
	 * @param b
	 *            bits in b determine which axes flip, 0 for no flip
	 * @return {@link SecondaryStructProjection} on the norm vectors
	 */
	public SecondaryStructProjection getAlphaNormProjection(byte b) {
		return alphaHelices.getNormProjection(b);
	}

	/**
	 * returns the norm projections for the beta strands using the bits in b to determine which axes to flip
	 * 
	 * @param b
	 *            bits in b determine which axes flip, 0 for no flip
	 * @return {@link SecondaryStructProjection} on the norm vectors
	 */
	public SecondaryStructProjection getBetaNormProjection(byte b) {
		return betaStrands.getNormProjection(b);
	}

	/**
	 * Prints a {@link SecondaryStructProjection}
	 * 
	 * @param p
	 *            {@link SecondaryStructProjection} to print
	 */
	public static void printProjection(SecondaryStructProjection p) {
		for (int i = 0; i < p.length(); i++) {
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f" + System.lineSeparator(), p.getStart(i).x, p.getStart(i).y,
					p.getEnd(i).x, p.getEnd(i).y);
		}
	}

	/**
	 * prints the coordinates of the vectors
	 * 
	 * @param t
	 *            {@link Tuple2} of int[] for the vectors
	 */
	public void printVectors(Tuple2<int[], int[]> t) {
		int[] s = t._1;
		int[] e = t._2;
		for (int i = 0; i < s.length; i++) {
			System.out.println("=" + pts[s[i]] + "\t=" + pts[e[i]]);// + "\t=Segment[A" + (i + 1) + ", B" + (i + 1) +
			// "]");
		}
	}

	/**
	 * prints the alpha projection onto the ind'th alpha helix
	 * 
	 * @param ind
	 *            index of which alpha helix to be the plane for projection
	 */
	public void printAlphaProjection(int ind) {
		SecondaryStructProjection alpha = getAlphaProjection(ind);
		for (int i = 0; i < alpha.length(); i++) {
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f" + System.lineSeparator(), alpha.getStart(i).x,
					alpha.getStart(i).y, alpha.getEnd(i).x, alpha.getEnd(i).y);
		}
	}

	/**
	 * prints the beta projection onto the ind'th beta strand
	 * 
	 * @param ind
	 *            index of which beta strand to be the plane for projection
	 */
	public void printBetaProjection(int ind) {
		SecondaryStructProjection beta = getBetaProjection(ind);
		for (int i = 0; i < beta.length(); i++) {
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f" + System.lineSeparator(), beta.getStart(i).x,
					beta.getStart(i).y, beta.getEnd(i).x, beta.getEnd(i).y, i + 1);
		}
	}

	/**
	 * shows the alpha helices
	 */
	public void printHelices() {
		Tuple2<int[], int[]> helices = alphaHelices.getFeatures();
		for (int i = 0; i < helices._1.length; i++) {
			System.out.printf("%d-%d, %.5f" + System.lineSeparator(), helices._1[i] + 1, helices._2[i] + 1,
					pts[helices._1[i]].distance(pts[helices._2[i]]));
		}
	}

	/**
	 * shows the beta strands
	 */
	public void printStrands() {
		Tuple2<int[], int[]> strand = betaStrands.getFeatures();
		for (int i = 0; i < strand._1.length; i++) {
			boolean nice = false;
			if (nice)
				System.out.printf("%d-%d|", strand._1[i] + 1, strand._2[i] + 1);
			else
				System.out.printf("%d-%d, %.5f :: %s" + System.lineSeparator(), strand._1[i] + 1, strand._2[i] + 1,
						pts[strand._1[i]].distance(pts[strand._2[i]]), pts[strand._1[i]]);
		}
	}

	/**
	 * shows the points and the x y z coordinates
	 */
	public void printPoints() {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < pts.length; i++) {
			s.append(i);
			s.append(":\t");
			s.append(pts[i]);
			s.append(System.lineSeparator());
		}
		System.out.println(s);
	}

	/**
	 * Returns a {@link SecondaryStructureSequenceFeature} showing the alpha helices and beta strands.
	 * 
	 * @return A {@link SecondaryStructureSequenceFeature} showing the alpha helices and beta strands.
	 */
	public SecondaryStructureSequenceFeature getSequenceFeature() {
		Tuple2<int[], int[]> a, b;
		a = SecondaryStructTools.alphaHelices(dists, pts, ALPHA_FILTER, smooth);
		b = SecondaryStructTools.betaStrands(dists, pts, BETA_FILTER, smooth);
		return new SecondaryStructureSequenceFeature(dists.length, a, b);
	}
}
