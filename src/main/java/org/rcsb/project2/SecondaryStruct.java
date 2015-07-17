package org.rcsb.project2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import javax.vecmath.Point3d;
import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;

import scala.Tuple2;

public class SecondaryStruct {
	private static final int NO_MATCH_PENALTY = 50;
	public static final double C = 3;
	public static final double C_SQ = C * C;

	private static final int NUM = 4;

	// Feature that identifies beta strands
	private static final Feature BETA_STRAND = new Feature() {
		private static final double b2 = 6.2; // distance of points 2 apart
		private static final double b3 = 9.5; // distance of points 3 apart
		private static final double b4 = 12.0; // distance of points 4 apart
		private static final double bDiffThreshold = 0.21;

		@Override
		public boolean match(int i, double d) {
			if (i < 2)
				return true;
			if (i == 2)
				return d > b2;
			if (i == 3)
				return d > b3;
			if (i == 4)
				return d > b4;
			return false;
		}

		@Override
		public boolean match(double[] d) {
			double diff = 0;
			diff += Math.max((b2 - d[1]), 0);
			diff += Math.max((b3 - d[2]), 0);
			diff += Math.max((b4 - d[3]), 0);
			return diff <= bDiffThreshold;
		}
	};

	// Feature that identifies Alpha Helices
	private static final Feature ALPHA_HELIX = new Feature() {
		private static final double a2l = 5.0;
		private static final double a2h = 6.0;
		private static final double a3l = 4.0;
		private static final double a3h = 6.0;
		private static final double a4l = 5.0;
		private static final double a4h = 7.0;

		@Override
		public boolean match(int i, double d) {
			if (i < 2)
				return true;
			if (i == 2)
				return a2l < d && d < a2h;
			if (i == 3)
				return a3l < d && d < a3h;
			if (i == 4)
				return a4l < d && d < a4h;
			return false;
		}

		@Override
		public boolean match(double[] d) {
			return a2l < d[1] && d[1] < a2h && a3l < d[2] && d[2] < a3h && a4l < d[3] && d[3] < a4h;
		}

	};

	private Point3d[] pts;
	private double[][] dists;
	private Tuple2<int[], int[]> alphaHelices = null;
	private Tuple2<int[], int[]> betaStrands = null;
	SecondaryStructProjection[] alphaProjections = null;
	SecondaryStructProjection[] betaProjections = null;

	public SecondaryStruct(Point3d[] pts) {
		this.pts = pts;
		dists = dists(pts);
		alphaHelices = alphaHelices(dists, 4);
		betaStrands = betaStrands(dists, 3);
		alphaProjections = new SecondaryStructProjection[getAlphaLength()];
		betaProjections = new SecondaryStructProjection[getBetaLength()];
	}

	/**
	 * The number of points in this protein chain.
	 * 
	 * @return The number of points in this protein chain.
	 */
	public int length() {
		return dists.length;
	}

	public int getAlphaLength() {
		return alphaHelices._1.length;
	}

	public int getBetaLength() {
		return betaStrands._1.length;
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

	public Tuple2<int[], int[]> getHelices() {
		return alphaHelices;
	}

	public Tuple2<int[], int[]> getStrands() {
		return betaStrands;
	}

	/**
	 * Return a {@link SecondaryStructProjection} for the projection of all the alpha helix vectors onto the ith alpha
	 * helix plane
	 * 
	 * @param i
	 * @return
	 */
	public SecondaryStructProjection getAlphaProjection(int i) {
		return alphaProjections[i] == null ? alphaProjections[i] = project(pts, getHelices(), i) : alphaProjections[i];
	}

	public SecondaryStructProjection getBetaProjection(int i) {
		return betaProjections[i] == null ? betaProjections[i] = project(pts, getStrands(), i) : betaProjections[i];
	}

	public void printAlphaProjection(int ind) {
		SecondaryStructProjection alpha = getAlphaProjection(ind);
		for (int i = 0; i < alpha.length(); i++) {
			System.out.printf("%.3f\t%.3f\t%.3f\t%.3f" + System.lineSeparator(), alpha.getStart(i).x,
					alpha.getStart(i).y, alpha.getEnd(i).x, alpha.getEnd(i).y);
		}
	}

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
		Tuple2<int[], int[]> helices = getHelices();
		for (int i = 0; i < helices._1.length; i++) {
			System.out.printf("%d-%d" + System.lineSeparator(), helices._1[i] + 1, helices._2[i] + 1);
		}
	}

	/**
	 * shows the beta strands
	 */
	public void printStrands() {
		Tuple2<int[], int[]> strand = getStrands();
		for (int i = 0; i < strand._1.length; i++) {
			System.out.printf("%d-%d" + System.lineSeparator(), strand._1[i] + 1, strand._2[i] + 1);
		}
	}

	/**
	 * shows the points and the x y z coordinates
	 */
	public void printPoints() {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < pts.length; i++) {
			s.append(i + 1);
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
		a = alphaHelices(dists, 4);
		b = betaStrands(dists, 3);
		return new SecondaryStructureSequenceFeature(dists.length, a, b);
	}

	/**
	 * Computes the distances for the given array of points.
	 * 
	 * @param pts
	 *            Array of {@link Point3d} representing the points of the protein chain.
	 * @return The distances for the given array of points.
	 */
	public static double[][] dists(Point3d[] pts) {
		double[][] out = new double[pts.length][NUM];
		for (int i = 0; i < NUM; i++) {
			for (int j = 0; j < i; j++) {
				if (pts[i] != null && pts[i - j - 1] != null)
					out[i][j] = pts[i].distance(pts[i - j - 1]);
			}
		}
		for (int i = NUM; i < pts.length; i++) {
			for (int j = 0; j < NUM; j++) {
				if (pts[i] != null && pts[i - j - 1] != null)
					out[i][j] = pts[i].distance(pts[i - j - 1]);
			}
		}
		return out;
	}

	/**
	 * Makes the distances easier to see
	 * 
	 * @param dists
	 *            Array of distances.
	 * @return String showing the distances.
	 */
	public static String distsToString(double[][] dists) {
		return distsToString(dists, 0);
	}

	/**
	 * Makes the distances easier to see
	 * 
	 * @param dists
	 *            Array of distances.
	 * @param offset
	 *            integer offset to start count from.
	 * @return String showing the distances.
	 */
	public static String distsToString(double[][] dists, int offset) {
		if (dists == null)
			return null;
		StringBuilder a;
		a = new StringBuilder();
		for (int i = 0; i < dists.length; i++)
			a.append((i + offset) + "\t");
		a.append(System.lineSeparator());
		for (int j = 0; j < NUM; j++) {
			for (int i = 0; i < dists.length; i++)
				a.append(String.format("%.3f", dists[i][j]) + "\t");
			a.append(System.lineSeparator());
		}
		return a.toString();
	}

	/**
	 * Main method for matching features
	 * 
	 * @param feat
	 *            Feature to be matched to.
	 * @param dists
	 *            Array of distances.
	 * @param filter
	 *            Any sequence of points that matches to the feature with less than this many points is filtered out.
	 * @return @{link Tuple2} of 2 arrays, 1 is starting indeces, 2 is ending indeces. Some may overlap.
	 */
	public static Tuple2<int[], int[]> match(Feature feat, double[][] dists, int filter) {
		List<Integer> s = new ArrayList<>();
		List<Integer> e = new ArrayList<>();
		int on = 1;
		out: for (int i = 0; i < dists.length; i++) {
			if (on == NUM) {
				// System.out.println("Start : " + i);
				// System.out.println("Real Start : " + (i - NUM));
				s.add(i - NUM);
				while (i < dists.length && feat.match(dists[i]))
					i++;
				// System.out.println(i + " breaks");
				e.add(i - 1);
				on--;
				i--;
				continue;
			}
			for (int j = 0; j < on; j++) {
				// System.out.println(i + ": " + on);
				if (!feat.match(j + 1, dists[i][j])) {
					// System.out.printf("[%d][%d] is bad" + System.lineSeparator(), i, j);
					// System.out.println("on is now " + (j + 1));
					on = j + 1;
					continue out;
				}
				// System.out.printf("[%d][%d] is ok" + System.lineSeparator(), i, j);
			}
			on++;
			// System.out.println("on reached " + on);
		}
		int size = Math.min(s.size(), e.size());
		if (s.size() != e.size())
			new Exception("s and e are not the same size: " + s.size() + ", " + e.size() + ". Will trim to match.")
					.printStackTrace();
		for (int i = 0; i < size; i++) {
			while (i < e.size() && i < s.size() && e.get(i) - s.get(i) <= filter) {
				s.remove(i);
				e.remove(i);
			}
		}
		size = Math.min(e.size(), s.size());
		int[] so = new int[size];
		int[] eo = new int[size];
		for (int i = 0; i < size; i++) {
			so[i] = s.get(i);
			eo[i] = e.get(i);
		}
		return new Tuple2<>(so, eo);
	}

	/**
	 * Returns a {@link Tuple2} of int[] showing where the alpha helices start and end.
	 * 
	 * @param dists
	 *            Array of distances.
	 * @param filter
	 *            Any sequence of points that matches to the feature with less than this many points is filtered out.
	 * @return @{link Tuple2} of 2 arrays, 1 is starting indeces, 2 is ending indeces. Some may overlap.
	 */
	public static Tuple2<int[], int[]> alphaHelices(double[][] dists, int filter) {
		return match(ALPHA_HELIX, dists, filter);
	}

	/**
	 * Returns a {@link Tuple2} of int[] showing where the beta strands start and end.
	 * 
	 * @param dists
	 *            Array of distances.
	 * @param filter
	 *            Any sequence of points that matches to the feature with less than this many points is filtered out.
	 * @return @{link Tuple2} of 2 arrays, 1 is starting indeces, 2 is ending indeces. Some may overlap.
	 */
	public static Tuple2<int[], int[]> betaStrands(double[][] dists, int filter) {
		return match(BETA_STRAND, dists, filter);
	}

	public static SecondaryStructProjection project(Point3d[] pts, Tuple2<int[], int[]> d, int ind) {
		int[] s = d._1;
		int[] e = d._2;
		if (s.length != e.length)
			throw new IllegalArgumentException("Lengths do not match");
		int N = s.length;
		if (N == 1)
			return new SecondaryStructProjection(new Vector2d[] { new Vector2d(0, 0) }, new Vector2d[] { new Vector2d(
					0, 0) });
		Point3d start = pts[s[ind]];
		Point3d end = pts[e[ind]];
		Vector3d v = new Vector3d(end);
		v.sub(start);
		Vector3d x, y;
		x = new Vector3d(pts[s[(ind + 1) % N]]);
		if (x.equals(end))
			x = new Vector3d(pts[s[(ind + 2) % N]]);
		x.sub(end);
		y = new Vector3d();
		y.cross(x, v);
		x.normalize();
		y.normalize();
		Vector2d[] vecs = new Vector2d[N];
		Vector2d[] vece = new Vector2d[N];
		for (int i = 0; i < N; i++) {
			Vector3d vs = new Vector3d(pts[s[i]]);
			vs.sub(end);
			vecs[i] = projectPlane(vs, v, x, y);

			Vector3d ve = new Vector3d(pts[e[i]]);
			ve.sub(end);
			vece[i] = projectPlane(ve, v, x, y);
		}
		return new SecondaryStructProjection(vecs, vece);
	}

	public static Vector2d projectPlane(Vector3d vec, Vector3d plane, Vector3d x, Vector3d y) {
		Vector3d vect = new Vector3d(vec);
		vect.sub(project(vect, plane));
		return new Vector2d(vect.dot(x), vect.dot(y));
	}

	public static Vector3d project(Vector3d vec, Vector3d onto) {
		Vector3d o = new Vector3d(onto);
		o.scale(vec.dot(onto) / onto.lengthSquared());
		return o;
	}

	public static double simil(Tuple2<Vector2d, Vector2d> v1, Tuple2<Vector2d, Vector2d> v2) {
		return simil(v1._1, v1._2, v2._1, v2._2);
	}

	public static double simil(Vector2d s1, Vector2d e1, Vector2d s2, Vector2d e2) {
		Vector2d ds = new Vector2d(s1);
		ds.sub(s2);
		Vector2d de = new Vector2d(e1);
		de.sub(e2);
		return ds.dot(ds) + de.dot(de) - de.dot(ds);
	}

	public static Tuple2<Double, int[]> rmsd(SecondaryStructProjection p1, SecondaryStructProjection p2) {
		double o1, o2;
		o1 = o2 = 0;
		int no1, no2;
		no1 = no2 = 0;

		int n1 = p1.length();
		int n2 = p2.length();
		int[] m1, m2;
		Tuple2<int[], int[]> t = align(p1, p2);
		m1 = t._1;
		m2 = t._2;
		// int[] m1 = align(p1, p2);
		// int[] m2 = align(p2, p1);
		// System.out.println(Arrays.toString(m1));
		// System.out.println(Arrays.toString(m2));
		for (int i = 0; i < m1.length; i++)
			if (m1[i] != -1) {
				double sim = simil(p1.get(i), p2.get(m1[i]));
				// System.out.println(i + ", " + m1[i] + ": " + sim);
				o1 += sim;
			}
			else {
				no1++;
				// System.out.println(i + ": " + NO_MATCH_PENALTY);
			}
		for (int i = 0; i < m2.length; i++)
			if (m2[i] != -1) {
				double sim = simil(p2.get(i), p1.get(m2[i]));
				// System.out.println(i + ", " + m2[i] + ": " + sim);
				o2 += sim;
			}
			else {
				no2++;
				// System.out.println(i + ": " + NO_MATCH_PENALTY);
			}
		// System.out.println("o1: " + o1);
		// System.out.println("o2: " + o2);
		o1 += no1 * NO_MATCH_PENALTY * (no1) / (double) n1;
		o2 += no2 * NO_MATCH_PENALTY * (no2) / (double) n2;
		// System.out.println("o1: " + o1);
		// System.out.println("o2: " + o2);
		o1 /= n1;
		o2 /= n2;
		o1 = Math.sqrt(o1);
		o2 = Math.sqrt(o2);
		// System.out.println("o1: " + o1);
		// System.out.println("o2: " + o2);
		return new Tuple2<>(Math.min(o1, o2), m1);
	}

	public static Tuple2<int[], int[]> align(SecondaryStructProjection p1, SecondaryStructProjection p2) {
		final double MIN_SIM = 3.0;
		int n1 = p1.length();
		int n2 = p2.length();
		int[] m1 = new int[n1];
		int[] m2 = new int[n2];
		// int yay = 0;
		for (int i = 0; i < n1; i++) {
			// System.out.print(i + " + ");
			Tuple2<Integer, Double> t = p2.getCloseTo(p1.get(i));
			m1[i] = t._1;
			if (t._2 != -1 && t._2 < MIN_SIM) {
				// System.out.println("yay " + ++yay);
				m2[m1[i]] = i;
			}
		}
		// System.out.println(Arrays.toString(m2));
		int test = 0;
		for (int i = 0; i < n2; i++) {
			if (m2[i] == 0) {
				// System.out.println("hi " + ++test);
				m2[i] = p1.getCloseTo(p2.get(i))._1;
			}
		}
		return new Tuple2<>(m1, m2);
	}

	/**
	 * Writes {@link Point3d} that goes with the specified pdbID array to a file to be read later.
	 * 
	 * @param pts
	 *            Array of {@link Point3d} to be written.
	 * @param name
	 *            The pdbID that corresponds to the points. This will be the filename.
	 */
	public static void write(Point3d[] pts, String name) {
		try (PrintWriter pw = new PrintWriter(new FileWriter(name + ".txt"))) {
			pw.println(pts.length);
			for (Point3d p : pts)
				if (p == null)
					pw.println();
				else
					pw.printf("%.3f %.3f %.3f" + System.lineSeparator(), p.x, p.y, p.z);
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Reads the Array of {@link Point3d} stored by using {@link SecondaryStruct#write }.
	 * 
	 * @param name
	 *            pdbID to read the points of.
	 * @return Array of {@link Point3d}.
	 */
	public static Point3d[] read(String name) {
		Point3d[] pts = null;
		try (Scanner scan = new Scanner(new File(name + ".txt"))) {
			pts = new Point3d[Integer.parseInt(scan.nextLine())];
			for (int i = 0; i < pts.length; i++) {
				String line = scan.nextLine();
				if (line.length() < 2)
					pts[i] = null;
				else {
					String[] spl = line.split(" ");
					pts[i] = new Point3d(Double.parseDouble(spl[0]), Double.parseDouble(spl[1]),
							Double.parseDouble(spl[2]));
				}
			}
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return pts;
	}

	/**
	 * Downloads the protein chain with given pdbID and returns an array of {@link Point3d} for the points in the chain.
	 * 
	 * @param pdbID
	 *            pdbID of the protein chain to download.
	 * @return Array of {@link Point3d} for the points in the chain.
	 * @throws IOException
	 * @throws StructureException
	 */
	public static Point3d[] pull(String pdbID) throws IOException, StructureException {
		Structure str;
		Atom[] atoms;
		str = StructureIO.getStructure(pdbID);
		atoms = StructureTools.getAtomCAArray(str);
		Point3d[] o = new Point3d[atoms.length];
		for (int i = 0; i < o.length; i++)
			o[i] = new Point3d(atoms[i].getX(), atoms[i].getY(), atoms[i].getZ());
		return o;
	}

	public static void align(SecondaryStruct s1, SecondaryStruct s2, String str1, String str2, boolean alpha, int i1,
			int i2) {
		System.out.println(str1 + " vs " + str2);
		System.out.println(alpha ? "Alpha Helices" : "Beta Strands");
		System.out.println(SecondaryStruct.rmsd(alpha ? s1.getAlphaProjection(i1) : s1.getBetaProjection(i1),
				alpha ? s2.getAlphaProjection(i2) : s2.getBetaProjection(i2))._1);
	}

	public static void align(SecondaryStruct s1, SecondaryStruct s2, String str1, String str2, boolean alpha) {
		final double MIN_PERFECT_RMSD = 3.0;
		System.out.println(str1 + " vs " + str2);
		System.out.println(alpha ? "Alpha Helices" : "Beta Strands");
		double minRMSD = Double.MAX_VALUE;
		int minInd1, minInd2;
		minInd1 = minInd2 = 0;
		int n1, n2;
		n1 = alpha ? s1.getAlphaLength() : s1.getBetaLength();
		n2 = alpha ? s2.getAlphaLength() : s2.getBetaLength();
		for (int i1 = 0; i1 < n1; i1++) {
			for (int i2 = 0; i2 < n2; i2++) {
				Tuple2<Double, int[]> t = SecondaryStruct.rmsd(
						alpha ? s1.getAlphaProjection(i1) : s1.getBetaProjection(i1), alpha ? s2.getAlphaProjection(i2)
								: s2.getBetaProjection(i2));
				double tes = t._1;

				if (tes < minRMSD) {
					minRMSD = tes;
					minInd1 = i1;
					minInd2 = i2;
				}
				// System.out.println("RMSD " + i1 + ", " + i2 + ": " + tes);
			}
		}
		System.out.println(minInd1 + "," + minInd2 + " min RMSD: " + minRMSD);
	}
}
