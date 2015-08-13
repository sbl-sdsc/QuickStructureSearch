package org.rcsb.project2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;
import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;

import scala.Tuple2;
import scala.Tuple3;

public class SecondaryStructTools {

	// private static final int MERGE = 5;
	private static final double ANGLE_DIFF = 2.4;
	private static final double MIN_REDUCE = 0.25;

	/**
	 * Calculates the similarity score between 2 protein chains, lower is more similar
	 * 
	 * @param p1
	 *            Array of {@link Point3d} for the first protein chain
	 * @param p2
	 *            Array of {@link Point3d} for the second protein chain
	 * @param smooth
	 *            whether the points were smoothed or not
	 * @return similarity score between 2 protein chains
	 */
	public static double calculateScore(Point3d[] p1, Point3d[] p2, boolean smooth) {
		SecondaryStruct s1 = new SecondaryStruct(p1, smooth);
		SecondaryStruct s2 = new SecondaryStruct(p2, smooth);
		return align(s1, s2);
	}

	public static double align(SecondaryStruct s1, SecondaryStruct s2) {
		// if (s1.getAlphaLength() < MERGE && s2.getAlphaLength() < MERGE || s1.getBetaLength() < MERGE
		// && s2.getBetaLength() < MERGE) {
		// System.out.println("MERGED");
		// return align(new MergeStruct(s1), new MergeStruct(s2));
		// }
		double a = align(s1.getAlpha(), s2.getAlpha());
		double b = align(s1.getBeta(), s2.getBeta());
		// System.out.println("TM: " + TmScorer.getFatCatTmScore(s1.getPoints(), s2.getPoints())[0]);
		// System.out.println("A: " + a);
		// System.out.println("B: " + b);
		// if (s1.getAlphaLength() == 1 || s2.getAlphaLength() == 1 || s1.getBetaLength() == 1 || s2.getBetaLength() ==
		// 1)
		// return -1;
		int al = Math.min(s1.getNumAlphaPoints(), s2.getNumAlphaPoints());
		int bl = Math.min(s1.getNumBetaPoints(), s2.getNumBetaPoints());
		// System.out.println("Alen: " + al);
		// System.out.println("Blen: " + bl);
		double c = (a * (al) + b * (bl)) / (al + bl);
		return c;
	}

	/**
	 * Aligns two {@link StructureCompare}s
	 * 
	 * @param f1
	 *            first {@link StructureCompare}
	 * @param f2
	 *            second {@link StructureCompare}
	 * @return A double for how similar the two are. Lower is more similar, higher is less similar
	 */
	public static double align(StructureCompare f1, StructureCompare f2) {
		final boolean print = false;
		final double MIN_PERFECT_RMSD = 2.0; // if some alignment scores below this, then stop trying other alignments.
		int[] perfectMatch = new int[0]; // Array storing the perfect matched pairs
		// if some alignment scores lower than MIN_PERFECT_RMSD then only match pairs of the perfect match
		double minRMSD = Double.MAX_VALUE;
		for (byte b = 0; b < 4; b++) {// 0-4 covers all orientations
			Tuple2<Double, int[]> t = rmsd(f1.getNormProjection((byte) (b % 4)),// probably can replace with just b
					f2.getNormProjection((byte) (b >> 2))); // probably can replace this with just 0
			if (print) {
				System.out.println(b + " : " + t._1);
				System.out.println(Arrays.toString(t._2));
			}
			minRMSD = Math.min(minRMSD, t._1);
			if (minRMSD < MIN_PERFECT_RMSD) {
				perfectMatch = t._2;
				break;
			}
		}
		if (print) {
			System.out.println("NormP: " + minRMSD);
			System.out.println("Perfect: " + Arrays.toString(perfectMatch));
		}
		int n1, n2;
		n1 = f1.length(); // number of features
		n2 = f2.length();
		// System.out.println("n1 : " + n1 + ", n2: " + n2);
		if (n1 == 0 || n2 == 0)
			return 10;
		// loop through all pairs
		outer: for (int i1 = 0; i1 < n1; i1++)
			for (int i2 = 0; i2 < n2; i2++) {
				if (minRMSD < MIN_PERFECT_RMSD) { // if the minRMSD is low enough, only check perfect match
					for (int i = i1; i < perfectMatch.length; i++)
						if (perfectMatch[i] != -1) {
							double test = rmsd(f1.getProjection(i), f2.getProjection(perfectMatch[i]))._1;
							if (print)
								System.out.println("Perfect Match minRMSD " + i + ", " + perfectMatch[i] + ": " + test);
							if (test < minRMSD)
								minRMSD = test;
						}
					break outer;
				}
				// minRMSD not low enough, check the projection of f1 and f2 onto i1 and i2 respectively
				Tuple2<Double, int[]> t = rmsd(f1.getProjection(i1), f2.getProjection(i2));
				double tes = t._1;
				if (tes < minRMSD) {
					minRMSD = tes;
					perfectMatch = t._2;
					if (print) {
						System.out.println("ALN: " + Arrays.toString(t._2));
						System.out.println("minRMSD " + i1 + ", " + i2 + ": " + tes);
					}
				}
				if (print)
					System.out.println("RMSD " + i1 + ", " + i2 + ": " + tes);
			}
		return minRMSD;
	}

	/**
	 * Matches {@link SecondaryStruct} projections and prints
	 * 
	 * @param s1
	 *            First {@link SecondaryStruct}
	 * @param s2
	 *            Second {@link SecondaryStruct}
	 * @param str1
	 *            Name of first protein chain
	 * @param str2
	 *            Name of second protein chain
	 * @param alpha
	 *            if true matches alpha helices, else matches beta strands
	 * @param i1
	 *            index to project the first {@link SecondaryStruct}
	 * @param i2
	 *            index to project the second {@link SecondaryStruct}
	 */
	public static void align(SecondaryStruct s1, SecondaryStruct s2, String str1, String str2, boolean alpha, int i1,
			int i2) {
		System.out.println(str1 + " vs " + str2);
		System.out.println(alpha ? "Alpha Helices" : "Beta Strands");
		System.out.println(SecondaryStructTools.rmsd(alpha ? s1.getAlphaProjection(i1) : s1.getBetaProjection(i1),
				alpha ? s2.getAlphaProjection(i2) : s2.getBetaProjection(i2))._1);
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

	/**
	 * Reads the Array of {@link Point3d} stored by using {@link SecondaryStructTools#write }.
	 * 
	 * @param name
	 *            pdbID to read the points of.
	 * @return Array of {@link Point3d}.
	 */
	public static Point3d[] read(String name) {
		Point3d[] pts = null;
		try (BufferedReader reader = new BufferedReader(new FileReader("data/" + name + ".txt"))) {
			pts = new Point3d[Integer.parseInt(reader.readLine())];
			for (int i = 0; i < pts.length; i++) {
				String line = reader.readLine();
				if (line.length() < 2)
					pts[i] = null;
				else {
					String[] spl = line.split(" ");
					pts[i] = new Point3d(Double.parseDouble(spl[0]), Double.parseDouble(spl[1]),
							Double.parseDouble(spl[2]));
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		return pts;
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
		try (PrintWriter pw = new PrintWriter(new FileWriter("data/" + name + ".txt"))) {
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
	 * If file doesnt exist, downloads and write to file, otherwise reads file
	 * 
	 * @param chain
	 *            Name of protein chain
	 * @return Point3d[] for the points in the protein chain
	 */
	public static Point3d[] obtain(String chain) {
		File f = new File("data/" + chain + ".txt");
		if (f.exists())
			return read(chain);
		Point3d[] pts = null;
		try {
			pts = SecondaryStructTools.pull(chain);
		}
		catch (IOException | StructureException e) {
			e.printStackTrace();
		}
		SecondaryStructTools.write(pts, chain);
		return pts;
	}

	/**
	 * Matches two {@link SecondaryStructProjection}s
	 * 
	 * @param p1
	 *            First {@link SecondaryStructProjection}
	 * @param p2
	 *            Second {@link SecondaryStructProjection}
	 * @return 3 int[] for matches in first to second, matches in second to first, and perfect matches
	 */
	public static Tuple3<int[], int[], int[]> align(SecondaryStructProjection p1, SecondaryStructProjection p2) {
		final double MIN_SIM = 3.0;
		int n1 = p1.length();
		int n2 = p2.length();
		int[] m1 = new int[n1];
		int[] m2 = new int[n2];
		int[] perf = new int[n1];
		for (int i = 0; i < n1; i++)
			m1[i] = perf[i] = -1;
		for (int i = 0; i < n2; i++)
			m2[i] = -1;
		// int yay = 0;
		for (int i = 0; i < n1; i++) {
			// System.out.print(i + " + ");
			Tuple2<Integer, Double> t = p2.getCloseTo(p1.get(i));
			if (t._1 != -1 && p1.getCloseTo(p2.get(t._1))._1 == i) {
				m1[i] = t._1;
				m2[m1[i]] = i;
			}
			if (t._2 != -1 && t._2 < MIN_SIM) {
				// System.out.println("yay " + ++yay);
				perf[i] = m1[i];
			}
		}
		// System.out.println(Arrays.toString(m2));
		// int test = 0;
		// for (int i = 0; i < n2; i++) {
		// if (m2[i] == 0) {
		// // System.out.println("hi " + ++test);
		// m2[i] = p1.getCloseTo(p2.get(i))._1;
		// }
		// }
		return new Tuple3<>(m1, m2, perf);
	}

	/**
	 * Calculates the score between two {@link SecondaryStructProjection}s
	 * 
	 * @param p1
	 *            First {@link SecondaryStructProjection}
	 * @param p2
	 *            Second {@link SecondaryStructProjection}
	 * @return Score and int[] for perfect match
	 */
	public static Tuple2<Double, int[]> rmsd(SecondaryStructProjection p1, SecondaryStructProjection p2) {
		boolean print = false;
		double o1, o2;
		o1 = o2 = 0;
		int no1, no2;
		no1 = no2 = 0;

		int n1 = p1.length();
		int n2 = p2.length();
		if (n1 == 0 || n2 == 0)
			return new Tuple2<>(-1.0, null);
		int[] m1, m2, perf;
		Tuple3<int[], int[], int[]> t = align(p1, p2); // find matching pairs
		m1 = t._1();
		m2 = t._2();
		perf = t._3();
		if (print) {
			System.out.println(Arrays.toString(m1));
			System.out.println(Arrays.toString(m2));
		}
		for (int i = 0; i < m1.length; i++)
			if (m1[i] != -1) {
				double sim = SecondaryStructTools.simil(p1.get(i), p2.get(m1[i])); // calculate score1
				if (print)
					System.out.println(i + ", " + m1[i] + ": " + sim);
				o1 += sim;
			}
			else {
				no1++;
				if (print)
					System.out.println(i + ": " + NO_MATCH_PENALTY);
			}
		for (int i = 0; i < m2.length; i++)
			if (m2[i] != -1) {
				double sim = SecondaryStructTools.simil(p2.get(i), p1.get(m2[i])); // calculate score2
				if (print)
					System.out.println(i + ", " + m2[i] + ": " + sim);
				o2 += sim;
			}
			else {
				no2++;
				if (print)
					System.out.println(i + ": " + NO_MATCH_PENALTY);
			}
		// combine scores
		if (print) {
			System.out.println("o1: " + o1);
			System.out.println("o2: " + o2);
		}
		o1 += no1 * SecondaryStructTools.NO_MATCH_PENALTY;
		o2 += no2 * SecondaryStructTools.NO_MATCH_PENALTY;
		if (print) {
			System.out.println("o1: " + o1);
			System.out.println("o2: " + o2);
		}
		o1 /= n1;
		o2 /= n2;
		if (print) {
			System.out.println("o1: " + o1);
			System.out.println("o2: " + o2);
		}
		o1 = Math.sqrt(o1);
		o2 = Math.sqrt(o2);
		if (print) {
			System.out.println("o1: " + o1);
			System.out.println("o2: " + o2);
		}
		o1 *= MIN_REDUCE + (1 - MIN_REDUCE) * no1 / n1;
		o2 *= MIN_REDUCE + (1 - MIN_REDUCE) * no2 / n2;
		if (print) {
			System.out.println("o1: " + o1);
			System.out.println("o2: " + o2);
			System.out.println();
		}
		double o = (o1 * n1 + o2 * n2) / (n1 + n2);
		return new Tuple2<>(o, perf);
	}

	/**
	 * Similarity of two projected vectors
	 * 
	 * @param s1
	 *            Start of first vector
	 * @param e1
	 *            End of first vector
	 * @param s2
	 *            Start of second vector
	 * @param e2
	 *            End of second vector
	 * @return double for how similar the two are, higher is less similar, lower is more similar
	 */
	public static double simil(Vector2d s1, Vector2d e1, Vector2d s2, Vector2d e2) {
		Vector2d ds = new Vector2d(s1);
		ds.sub(s2);
		Vector2d de = new Vector2d(e1);
		de.sub(e2);
		return ds.dot(ds) + de.dot(de) - de.dot(ds);
	}

	/**
	 * Similarity of two projected vectors
	 * 
	 * @param v1
	 *            Pair of {@link Vector2d} for the first vector
	 * @param v2
	 *            Pair of {@link Vector2d} for the second vector
	 * @return double for how similar the two are, higher is less similar, lower is more similar
	 */
	public static double simil(Tuple2<Vector2d, Vector2d> v1, Tuple2<Vector2d, Vector2d> v2) {
		return simil(v1._1, v1._2, v2._1, v2._2);
	}

	/**
	 * Projects a vector onto another vector
	 * 
	 * @param vec
	 *            Vector to project
	 * @param onto
	 *            Vector to project onto
	 * @return new Vector of the projection
	 */
	public static Vector3d project(Vector3d vec, Vector3d onto) {
		Vector3d o = new Vector3d(onto);
		o.scale(vec.dot(onto) / onto.lengthSquared());
		return o;
	}

	/**
	 * Projects a vector on the given plane with the given x and y directions.
	 * 
	 * @param vec
	 *            Vector to project
	 * @param plane
	 *            Plane to project onto
	 * @param x
	 *            direction of x
	 * @param y
	 *            direction of y
	 * @return 2d Vector that is the projection onto the plane
	 */
	public static Vector2d projectPlane(Vector3d vec, Vector3d plane, Vector3d x, Vector3d y) {
		Vector3d vect = new Vector3d(vec);
		vect.sub(project(vect, plane));
		return new Vector2d(vect.dot(x), vect.dot(y));
	}

	/**
	 * Projects vectors onto the plane defined by the ind'th vector. Used for projections onto the ind'th vector
	 * 
	 * @param pts
	 *            Array for all the points
	 * @param vectors
	 *            Vectors to project
	 * @param ind
	 *            Which vector to define the plane to project onto
	 * @return {@link SecondaryStructProjection} for the projection
	 */
	public static SecondaryStructProjection project(Point3d[] pts, Tuple2<int[], int[]> vectors, int ind) {
		int[] s = vectors._1;
		int[] e = vectors._2;
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

	/**
	 * Projects features onto a planeVec, using a given x direction. Used for norm projections
	 * 
	 * @param pts
	 *            Array of Point3d for the original points
	 * @param feats
	 *            Vectors to project
	 * @param planeVec
	 *            Plane to project onto
	 * @param end
	 *            Point in 3d space used for offset
	 * @param x
	 *            Direction of the x-axis
	 * @return {@link SecondaryStructProjection} for the projection
	 */
	public static SecondaryStructProjection projectOnto(Point3d[] pts, Tuple2<int[], int[]> feats, Vector3d planeVec,
			Point3d end, Vector3d x) {
		if (x.equals(end))
			throw new IllegalArgumentException("X is same as end");
		int[] s = feats._1;
		int[] e = feats._2;
		Vector3d xN = new Vector3d(x);
		Vector3d yN = new Vector3d();
		yN.cross(x, planeVec);
		xN.scale(1 / xN.length());
		yN.scale(1 / xN.length());
		// xN.normalize();
		// yN.normalize();
		int N = feats._1.length;
		Vector2d[] vecs = new Vector2d[N];
		Vector2d[] vece = new Vector2d[N];
		for (int i = 0; i < N; i++) {
			Tuple2<Vector3d, Vector3d> t = getFeatureStartEnd(s[i], e[i], pts);
			Vector3d vs = t._1;
			Vector3d ve = t._2;
			// Vector3d vs = new Vector3d(pts[s[i]]);
			vs.sub(end);
			vecs[i] = projectPlane(vs, planeVec, xN, yN);

			// Vector3d ve = new Vector3d(pts[e[i]]);
			ve.sub(end);
			vece[i] = projectPlane(ve, planeVec, xN, yN);
		}
		return new SecondaryStructProjection(vecs, vece);
	}

	/**
	 * Gets an approximate start and end point for a secondary structure. This is used to get a better approximation of
	 * the direction of the secondary structure instead of going direction from start to end
	 * 
	 * @param start
	 *            index of start point
	 * @param end
	 *            index of end point
	 * @param pts
	 *            Array of Points
	 * @return Two {@link Vector3d}s for the start and end of the feature
	 */
	public static Tuple2<Vector3d, Vector3d> getFeatureStartEnd(int start, int end, Point3d[] pts) {
		int mid = (start + end) / 2;
		Vector3d s = new Vector3d();
		Vector3d e = new Vector3d();
		for (int i = start; i <= mid; i++)
			if (pts[i] != null)
				s.add(pts[i]);
		for (int i = mid; i <= end; i++)
			if (pts[i] != null)
				e.add(pts[i]);
		s.scale(1 / (double) (mid - start + 1));
		e.scale(1 / (double) (end - mid + 1));
		Vector3d eDiff = new Vector3d(e);
		eDiff.sub(s);
		Vector3d relE;
		int i = 0;
		while (pts[end - i] == null)
			i++;
		relE = new Vector3d(pts[end - i]);
		relE.sub(s);
		Vector3d relS = new Vector3d(pts[start]);
		relS.sub(s);
		relS = project(relS, eDiff);
		relE = project(relE, eDiff);
		relS.add(s);
		relE.add(s);
		// System.out.println(pts[start]);
		// System.out.println(pts[end]);
		// System.out.println();
		// System.out.println(s);
		// System.out.println(e);
		// System.out.println();
		// System.out.println(relS);
		// System.out.println(relE);
		// System.out.println();
		// System.out.println();
		// System.out.println();
		return new Tuple2<>(relS, relE);
	}

	/**
	 * Computes the distances for the given array of points.
	 * 
	 * @param pts
	 *            Array of {@link Point3d} representing the points of the protein chain.
	 * @return The distances for the given array of points.
	 */
	public static double[][] dists(Point3d[] pts) {
		double[][] out = new double[pts.length][SecondaryStruct.NUM];
		for (int i = 0; i < SecondaryStruct.NUM; i++)
			for (int j = 0; j < i; j++)
				if (pts[i] != null && pts[i - j - 1] != null)
					out[i][j] = pts[i].distance(pts[i - j - 1]);
		for (int i = SecondaryStruct.NUM; i < pts.length; i++)
			for (int j = 0; j < SecondaryStruct.NUM; j++)
				if (pts[i] != null && pts[i - j - 1] != null)
					out[i][j] = pts[i].distance(pts[i - j - 1]);
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
		return SecondaryStructTools.distsToString(dists, 0);
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
		for (int j = 0; j < SecondaryStruct.NUM; j++) {
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
	 * @param pts
	 *            Array of points
	 * @param filter
	 *            Any sequence of points that matches to the feature with less than this many points is filtered out.
	 * @return @{link Tuple2} of 2 arrays, 1 is starting indeces, 2 is ending indeces. Some may overlap.
	 */
	public static Tuple2<int[], int[]> match(Feature feat, double[][] dists, Point3d[] pts, int filter) {
		List<Integer> s = new ArrayList<>();
		List<Integer> e = new ArrayList<>();
		int on = 1;
		out: for (int i = 0; i < dists.length; i++) {
			if (on == SecondaryStruct.NUM) {
				// System.out.println("Start : " + i);
				// System.out.println("Real Start : " + (i - NUM));
				s.add(i - SecondaryStruct.NUM);
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
		// System.out.println("S: " + s);
		// System.out.println("E: " + e);
		for (int i = 0; i < size; i++) {
			while (i < e.size() && i < s.size() && e.get(i) - s.get(i) <= filter) {
				s.remove(i);
				e.remove(i);
			}
		}
		int maxSeparate = 2;
		for (int i = 1; i < s.size(); i++) {
			while (i < s.size() && s.get(i) - e.get(i - 1) <= maxSeparate) {
				Vector3d v1 = new Vector3d(pts[s.get(i)]);
				v1.sub(pts[e.get(i)]);
				Vector3d v2 = new Vector3d(pts[e.get(i - 1)]);
				v2.sub(pts[s.get(i - 1)]);
				double diff = v1.angle(v2);
				// System.out.printf("(%d,%d) with (%d,%d) :: %.5f\t\t%.5f" + System.lineSeparator(), s.get(i - 1),
				// e.get(i - 1), s.get(i), e.get(i), diff, pts[s.get(i - 1)].distance(pts[e.get(i - 1)]));
				if (diff > ANGLE_DIFF) {
					// System.out.println("remove : " + i);
					e.remove(i - 1);
					s.remove(i);
				}
				else
					break;
			}
		}
		size = s.size();
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
	 * @param pts
	 *            Array of points.
	 * @param filter
	 *            Any sequence of points that matches to the feature with less than this many points is filtered out.
	 * @param smooth
	 * @return @{link Tuple2} of 2 arrays, 1 is starting indeces, 2 is ending indeces. Some may overlap.
	 */
	public static Tuple2<int[], int[]> alphaHelices(double[][] dists, Point3d[] pts, int filter, boolean smooth) {
		return match(smooth ? SMOOTH_ALPHA_HELIX : ALPHA_HELIX, dists, pts, filter);
	}

	/**
	 * Returns a {@link Tuple2} of int[] showing where the beta strands start and end.
	 * 
	 * @param dists
	 *            Array of distances.
	 * @param pts
	 *            Array of points.
	 * @param filter
	 *            Any sequence of points that matches to the feature with less than this many points is filtered out.
	 * @param smooth
	 * @return @{link Tuple2} of 2 arrays, 1 is starting indeces, 2 is ending indeces. Some may overlap.
	 */
	public static Tuple2<int[], int[]> betaStrands(double[][] dists, Point3d[] pts, int filter, boolean smooth) {
		return match(smooth ? SMOOTH_BETA_STRAND : BETA_STRAND, dists, pts, filter);
	}

	static final int NO_MATCH_PENALTY = 200;

	// Feature that identifies Alpha Helices
	static final Feature ALPHA_HELIX = new Feature() {
		private static final double a2l = 4.7;
		private static final double a2h = 6.21;
		private static final double a3l = 4.0;
		private static final double a3h = 6.21;
		private static final double a4l = 5.0;
		private static final double a4h = 7.3;
		private static final double bDiffThreshold = 0.21;

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
			double er = 0;
			if (d[1] > a2h)
				er += d[1] - a2h;
			else if (d[1] < a2l)
				er += a2l - d[1];
			if (d[2] > a3h)
				er += d[2] - a3h;
			else if (d[2] < a3l)
				er += a3l - d[2];
			if (d[3] > a4h)
				er += d[3] - a4h;
			else if (d[3] < a4l)
				er += a4l - d[3];
			return er < bDiffThreshold;
		}

	};

	// Feature that identifies Beta Strands
	static final Feature BETA_STRAND = new Feature() {
		private static final double b2 = 6.2; // distance of points 2 apart
		private static final double b3 = 9.5; // distance of points 3 apart
		private static final double b4 = 12.2; // distance of points 4 apart
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

	// Feature that identifies smoothed Alpha Helices
	static final Feature SMOOTH_ALPHA_HELIX = new Feature() {
		private static final double a2l = 0;
		private static final double a2h = 3.5;
		private static final double a3l = 0;
		private static final double a3h = 5.0;
		private static final double a4l = 0;
		private static final double a4h = 6.5;
		private static final double bDiffThreshold = 0.0;

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
			double er = 0;
			if (d[1] > a2h)
				er += d[1] - a2h;
			else if (d[1] < a2l)
				er += a2l - d[1];
			if (d[2] > a3h)
				er += d[2] - a3h;
			else if (d[2] < a3l)
				er += a3l - d[2];
			if (d[3] > a4h)
				er += d[3] - a4h;
			else if (d[3] < a4l)
				er += a4l - d[3];
			return er <= bDiffThreshold;
		}

	};

	// Feature that identifies smoothed Beta Strands
	static final Feature SMOOTH_BETA_STRAND = new Feature() {
		private static final double b2 = 6; // distance of points 2 apart
		private static final double b3 = 9; // distance of points 3 apart
		private static final double b4 = 12; // distance of points 4 apart
		private static final double bDiffThreshold = 0.0;

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
}
