package org.rcsb.project2;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function2;
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
import org.jblas.util.Random;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.SeqToChainMapper;

import scala.Tuple2;

/**
 * FingerprintMapper gives a visual representation of the alignment of two protein chains.
 * 
 * @author Kevin Wu
 */
public class FingerprintMapper_KevinWu {
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	private static final char MATCH_CHAR = '|';
	private static final char MATCH_FRAG_CHAR = '$';
	private static final String CA_NAME = "CA";
	private static final String GROUP_NAME = "GLU";
	private static final String separator = "\t";

	public static void main(String[] args) {
		String path = args[0];

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]").setAppName(
				FingerprintMapper_KevinWu.class.getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);

		long start = System.nanoTime();

		// read sequence file and map sequence length to an RDD
		List<Tuple2<String, Point3d[]>> list = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS * NUM_TASKS_PER_THREAD)
				.sample(false, 0.006)// sample
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate array> pairs
				.filter(new GapFilter(0, 0)) // filter chains with zero gap length and zero gaps
				// .filter(new org.apache.spark.api.java.function.Function<Tuple2<String, Point3d[]>, Boolean>() {
				// @Override
				// public Boolean call(Tuple2<String, Point3d[]> v1) throws Exception {
				// return v1._1().toUpperCase().equals("1SRR.A") || v1._1().toUpperCase().equals("1T37.A");
				// }
				// })// filter test
				.collect();
		List<Integer> l = new ArrayList<>();
		for (int i = 0; i < list.size() - 1; i++) {
			Point3d[] chain1 = list.get(i)._2;
			Point3d[] chain2 = list.get(i + 1)._2;
			if (chain1.length == chain2.length)
				l.add(i);
		}
		l.parallelStream().forEach(new Consumer<Integer>() {
			@Override
			public void accept(Integer i) {
				Point3d[] chain1 = list.get(i)._2;
				Point3d[] chain2 = list.get(i + 1)._2;
				String chainId1 = list.get(i)._1;
				String chainId2 = list.get(i + 1)._1;
				int[] f1 = new int[chain1.length];
				int[] f2 = new int[chain1.length];
				int m = 5;
				for (int k = 0; k < f1.length; k++) {
					f1[k] = Random.nextInt(m);
					f2[k] = Random.nextInt(m);
				}
				StringBuilder a1 = new StringBuilder(align(chain1, chain2, f1, f2,
						new Function2<Integer, Integer, Double>() {
							@Override
							public Double call(Integer v1, Integer v2) throws Exception {
								return 1 / (1d + (v1 ^ v2));
							}
						}));
				a1.append(System.lineSeparator());
				a1.append(System.lineSeparator());
				a1.append(align(chain1, chain2));
				a1.append(System.lineSeparator());
				a1.append(chainId1 + ", " + chainId2);
				System.out.println(a1);
			}
		});
		sc.close();
		System.out.println("Time: " + (System.nanoTime() - start) / 1E9 + " sec.");
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
	 *            Array of Point3d representing the first protein chain
	 * @param points2
	 *            Array of Point3d representing the second protein chain
	 * @return Optimal alignment, given by AFPChain's getOptAln()
	 */
	public static int[][][] getOpt(Point3d[] points1, Point3d[] points2) {
		Atom[] ca1 = getCAAtoms(points1);
		Atom[] ca2 = getCAAtoms(points2);

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
		return afp.getOptAln();
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
	 * @param p1
	 *            Array of Point3d representing the first protein chain
	 * @param p2
	 *            Array of Point3d representing the second protein chain
	 * @param f1
	 *            Array of T representing the first fragment chain
	 * @param f2
	 *            Array of T representing the second fragment chain
	 * @param sim
	 *            Function for the similarity between 2 fragments
	 * @return String showing the alignment
	 */
	public static <T> String align(Point3d[] p1, Point3d[] p2, T[] f1, T[] f2, Function2<T, T, Double> sim) {
		return align(getOpt(p1, p2), f1, f2, sim);
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
	 * @param optAln
	 *            Optimal alignment, given by AFPChain's getOptAln()
	 * @param f1
	 *            Array of T representing the first fragment chain
	 * @param f2
	 *            Array of T representing the second fragment chain
	 * @param sim
	 *            Function for the similarity between 2 fragments
	 * @return String showing the alignment
	 */
	public static <T> String align(int[][][] optAln, T[] f1, T[] f2, Function2<T, T, Double> sim) {
		class Comp implements FragmentComparable<Comp> {
			T t;

			Comp(T t) {
				this.t = t;
			}

			@Override
			public double compare(Comp c) {
				try {
					return sim.call(t, c.t).doubleValue();
				}
				catch (Exception e) {
					e.printStackTrace();
				}
				return 0;
			}

			@Override
			public String toString() {
				return t.toString();
			}

			@Override
			public boolean equals(Object o) {
				return t.equals(o);
			}
		}
		Comp[] fp1 = new Comp[f1.length];
		Comp[] fp2 = new Comp[f1.length];
		for (int i = 0; i < f1.length; i++) {
			fp1[i] = new Comp(f1[i]);
			fp2[i] = new Comp(f2[i]);
		}
		return align(optAln, fp1, fp2);
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
	 * @param p1
	 *            Array of Point3d representing the first protein chain
	 * @param p2
	 *            Array of Point3d representing the second protein chain
	 * @param f1
	 *            Array of T representing the first fragment chain
	 * @param f2
	 *            Array of T representing the second fragment chain
	 * @return String showing the alignment
	 */
	public static <T extends FragmentComparable<T>> String align(Point3d[] p1, Point3d[] p2, T[] f1, T[] f2) {
		return align(getOpt(p1, p2), f1, f2);
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
	 * @param optAln
	 *            Optimal alignment, given by AFPChain's getOptAln()
	 * @param f1
	 *            Array of T representing the first fragment chain
	 * @param f2
	 *            Array of T representing the second fragment chain
	 * @return String showing the alignment
	 */
	public static <T extends FragmentComparable<T>> String align(int[][][] optAln, T[] f1, T[] f2) {
		if (f1.length != optAln[0][0].length || f1.length != f2.length) {
			throw new IllegalArgumentException("Fragment lengths do not match: " + optAln[0][0].length + ", f1: "
					+ f1.length + ", f2: " + f2.length);
		}
		int[] i1 = optAln[0][0];
		int[] i2 = optAln[0][1];
		StringBuilder sf1 = new StringBuilder();
		StringBuilder si1 = new StringBuilder();
		StringBuilder mid = new StringBuilder();
		StringBuilder si2 = new StringBuilder();
		StringBuilder sf2 = new StringBuilder();
		for (int i = 0; i < f1.length; i++) {
			if (!(i1[i] == 0 && i2[i] == 0) || i == 0) {
				mid.append(MATCH_CHAR);
				si1.append(i1[i]);
				si2.append(i2[i]);
			}
			else {
				si1.append(String.format("(%d)", i + i1[0]));
				si2.append(String.format("(%d)", i + i2[0]));
			}
			if (f1[i].equals(f2[i]))
				mid.append(MATCH_FRAG_CHAR);
			else
				mid.append(String.format("%.2f", f1[i].compare(f2[i])));
			sf1.append(f1[i]);
			sf2.append(f2[i]);
			mid.append(separator);
			si1.append(separator);
			si2.append(separator);
			sf1.append(separator);
			sf2.append(separator);
		}
		sf1.append(System.lineSeparator());
		sf1.append(si1);
		sf1.append(System.lineSeparator());
		sf1.append(mid);
		sf1.append(System.lineSeparator());
		sf1.append(si2);
		sf1.append(System.lineSeparator());
		sf1.append(sf2);
		return sf1.toString();
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
	 * @param p1
	 *            Array of Point3d representing the first protein chain
	 * @param p2
	 *            Array of Point3d representing the second protein chain
	 * @param f1
	 *            Array of int representing the first fragment chain
	 * @param f2
	 *            Array of int representing the second fragment chain
	 * @param sim
	 *            Function for the similarity between 2 fragments
	 * @return String showing the alignment
	 */
	public static String align(Point3d[] p1, Point3d[] p2, int[] f1, int[] f2, Function2<Integer, Integer, Double> sim) {
		return align(getOpt(p1, p2), f1, f2, sim);
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
	 * @param optAln
	 *            Optimal alignment, given by AFPChain's getOptAln()
	 * @param f1
	 *            Array of int representing the first fragment chain
	 * @param f2
	 *            Array of int representing the first fragment chain
	 * @param sim
	 *            Function for the similarity between 2 fragments
	 * @return String showing the alignment
	 */
	public static String align(int[][][] optAln, int[] f1, int[] f2, Function2<Integer, Integer, Double> sim) {
		class IntComp implements FragmentComparable<IntComp> {
			private int i;

			IntComp(int i) {
				this.i = i;
			}

			@Override
			public double compare(IntComp t) {
				try {
					return sim.call(i, t.i).doubleValue();
				}
				catch (Exception e) {
					e.printStackTrace();
				}
				return 0;
			}

			@Override
			public String toString() {
				return Integer.toString(i);
			}

			@Override
			public boolean equals(Object o) {
				return o instanceof IntComp && i == ((IntComp) o).i;
			}
		}
		IntComp[] fp1 = new IntComp[f1.length];
		IntComp[] fp2 = new IntComp[f1.length];
		for (int i = 0; i < fp1.length; i++) {
			fp1[i] = new IntComp(f1[i]);
			fp2[i] = new IntComp(f2[i]);
		}
		return align(optAln, fp1, fp2);
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
	 *            Array of Point3d representing the first protein chain
	 * @param p2
	 *            Array of Point3d representing the second protein chain
	 * @return String showing the alignment
	 */
	public static String align(Point3d[] p1, Point3d[] p2) {
		return align(getOpt(p1, p2));
	}

	/**
	 * Gives a String that shows the alignment of two protein chains (no fragments). <br>
	 * <br>
	 * Row 1 has positions of points in chain 1<br>
	 * Row 2 has matching<br>
	 * {@value #MATCH_CHAR} means the points at the position match (based on fatcat) <br>
	 * Row 3 has positions of points in chain 2<br>
	 * 
	 * @param optAln
	 *            Optimal alignment, given by AFPChain's getOptAln()
	 * @return String showing the alignment
	 */
	public static String align(int[][][] optAln) {
		if (optAln.length != 1)
			new IllegalArgumentException("More than 1 perfect match").printStackTrace();
		int[] top = optAln[0][0];
		int[] bot = optAln[0][1];
		StringBuilder topI = new StringBuilder();
		StringBuilder match = new StringBuilder();
		StringBuilder botI = new StringBuilder();
		for (int i = 0; i < top.length; i++) {
			if (!(top[i] == 0 && bot[i] == 0) || i == 0) {
				match.append(MATCH_CHAR);
				topI.append(top[i]);
				botI.append(bot[i]);
			}
			else {
				topI.append(String.format("(%d)", i + top[0]));
				botI.append(String.format("(%d)", i + bot[0]));
			}
			topI.append(separator);
			match.append(separator);
			botI.append(separator);
		}
		topI.append(System.lineSeparator());
		topI.append(match);
		topI.append(System.lineSeparator());
		topI.append(botI);
		return topI.toString();
	}

	// copied from TMScorer.java
	private static Atom[] getCAAtoms(Point3d[] points) {
		int gaps = 0;
		for (Point3d p : points) {
			if (p == null) {
				gaps++;
			}
		}
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
}
