/*
 * $Id: SmithWatermanGotoh.java,v 1.10 2006/02/09 13:27:36 ahmed Exp $
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package org.rcsb.project3;

/**
 * This is the helper class for SmithWatermanGotohP3
 * An implementation of the Smith-Waterman algorithm with Gotoh's improvement
 * for biological local pairwise sequence alignment.
 * 
 * <strong>Recursive definition:</strong>
 * <ul>
 * <li>
 * <strong>Base conditions:</strong>
 * <ul>
 * <li><code>V(0, 0) = 0</code></li>
 * <li><code>V(i, 0) = E(i, 0) = W<sub>g</sub> + iW<sub>s</sub></code></li>
 * <li><code>V(0, j) = F(0, j) = W<sub>g</sub> + jW<sub>s</sub></code></li>
 * </ul>
 * </li>
 * <li>
 * <strong>Recurrence relation:</strong>
 * <ul>
 * <li><code>V(i, j) = max{E(i, j), F(i, j), G(i, j)}</code>, where:</li>
 * <li><code>G(i, j) = V(i - 1, j - 1) + similarity(S<sub>i</sub>, T<sub>j</sub>)</code></li>
 * <li><code>E(i, j) = max{E(i, j - 1) + W<sub>s</sub>, V(i, j - 1) + W<sub>g</sub> + W<sub>s</sub>}</code></li>
 * <li><code>F(i, j) = max{F(i - 1, j) + W<sub>s</sub>, V(i - 1, j) + W<sub>g</sub> + W<sub>s</sub>}</code></li>
 * </ul>
 * </ul> 
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public final class SmithWatermanGotoh {
	/**
	 * Hidden constructor
	 */
	private SmithWatermanGotoh() {
		super();
	}

	/**
	 * Aligns two sequences by Smith-Waterman (local)
	 * 
	 * @param s1
	 *            sequene #1 ({@link Sequence})
	 * @param s2
	 *            sequene #2 ({@link Sequence})
	 * @param matrix
	 *            scoring matrix ({@link Matrix})
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @return alignment object contains the two aligned sequences, the
	 *         alignment score and alignment statistics
	 * @see Sequence
	 * @see Matrix
	 */
	public static <T> Alignment<T> align(SequenceFeatureInterface<T> s1, SequenceFeatureInterface<T> s2,
			double o, double e) {

		int m = s1.length() + 1;
		int n = s2.length() + 1;

		byte[] pointers = new byte[m * n];

		// Initializes the boundaries of the traceback matrix to STOP.
		for (int i = 0, k = 0; i < m; i++, k += n) {
			pointers[k] = Directions.STOP;
		}
		for (int j = 1; j < n; j++) {
			pointers[j] = Directions.STOP;
		}

		short[] sizesOfVerticalGaps = new short[m * n];
		short[] sizesOfHorizontalGaps = new short[m * n];
		for (int i = 0, k = 0; i < m; i++, k += n) {
			for (int j = 0; j < n; j++) {
				sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
			}
		}

		Cell cell = SmithWatermanGotoh.construct(s1, s2, o, e,
				pointers, sizesOfVerticalGaps, sizesOfHorizontalGaps);
		Alignment<T> alignment = SmithWatermanGotoh.traceback(s1, s2,
				pointers, cell, sizesOfVerticalGaps, sizesOfHorizontalGaps);
		alignment.setOriginalSequence1(s1);
		alignment.setOriginalSequence2(s2);
		alignment.setOpen(o);
		alignment.setExtend(e);
		return alignment;
	}

	/**
	 * Constructs directions matrix for the traceback
	 * 
	 * @param s1
	 *            sequence #1
	 * @param s2
	 *            sequence #2
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @return The cell where the traceback starts.
	 */
	private static <T> Cell construct(SequenceFeatureInterface<T> s1, SequenceFeatureInterface<T> s2, double o, double e, byte[] pointers, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		int m = s1.length() + 1;
		int n = s2.length() + 1;

		double f; // score of alignment x1...xi to y1...yi if xi aligns to yi
		double[] g = new double[n]; // score if xi aligns to a gap after yi
		double h; // score if yi aligns to a gap after xi
		double[] v = new double[n]; // best score of alignment x1...xi to
		// y1...yi
		double vDiagonal;

		g[0] = Float.NEGATIVE_INFINITY;
		h = Float.NEGATIVE_INFINITY;
		v[0] = 0;

		for (int j = 1; j < n; j++) {
			g[j] = Float.NEGATIVE_INFINITY;
			v[j] = 0;
		}

		double similarityScore, g1, g2, h1, h2;

		Cell cell = new Cell();

		for (int i = 1, k = n; i < m; i++, k += n) {
			h = Float.NEGATIVE_INFINITY;
			vDiagonal = v[0];
			for (int j = 1, l = k + 1; j < n; j++, l++) {
				similarityScore = s1.similarity(s2, i-1, j-1);

				// Fill the matrices
				f = vDiagonal + similarityScore;

				g1 = g[j] - e;
				g2 = v[j] - o;
				if (g1 > g2) {
					g[j] = g1;
					sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
				} else {
					g[j] = g2;
				}

				h1 = h - e;
				h2 = v[j - 1] - o;
				if (h1 > h2) {
					h = h1;
					sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
				} else {
					h = h2;
				}

				vDiagonal = v[j];
				v[j] = maximum(f, g[j], h, 0);

				// Determine the traceback direction
				if (v[j] == 0) {
					pointers[l] = Directions.STOP;
				} else if (v[j] == f) {
					pointers[l] = Directions.DIAGONAL;
				} else if (v[j] == g[j]) {
					pointers[l] = Directions.UP;
				} else {
					pointers[l] = Directions.LEFT;
				}

				// Set the traceback start at the current cell i, j and score
				if (v[j] > cell.getScore()) {
					cell.set(i, j, v[j]);
				}
			}
		}
		return cell;
	}

	/**
	 * Returns the alignment of two sequences based on the passed array of
	 * pointers
	 * 
	 * @param s1
	 *            sequence #1
	 * @param s2
	 *            sequence #2
	 * @param m
	 *            scoring matrix
	 * @param cell
	 *            The cell where the traceback starts.
	 * @return {@link Alignment}with the two aligned sequences and alignment
	 *         score.
	 * @see Cell
	 * @see Alignment
	 */
	private static <T> Alignment<T> traceback(SequenceFeatureInterface<T> s1, SequenceFeatureInterface<T> s2,
			byte[] pointers, Cell cell, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		int n = s2.length() + 1;

		Alignment<T> alignment = new Alignment<T>();
		alignment.setScore(cell.getScore());

		int maxlen = s1.length() + s2.length(); // maximum length after the
		// aligned sequences

		Integer[] reversed1 = new Integer[maxlen]; // reversed sequence #1
		Integer[] reversed2 = new Integer[maxlen]; // reversed sequence #2

		int len1 = 0; // length of sequence #1 after alignment
		int len2 = 0; // length of sequence #2 after alignment

		int identity = 0; // count of identitcal pairs
		int similarity = 0; // count of similar pairs
		int gaps = 0; // count of gaps

		int c1, c2;

		int i = cell.getRow(); // traceback start row
		int j = cell.getCol(); // traceback start col
		int k = i * n;

		boolean stillGoing = true; // traceback flag: true -> continue & false
		// -> stop

		while (stillGoing) {
			try {
			switch (pointers[k + j]) {
			case Directions.UP:
				for (int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = --i;
					reversed2[len2++] = null;
					k -= n;
					gaps++;
				}
				break;
			case Directions.DIAGONAL:
				c1 = --i;
				c2 = --j;
				k -= n;
				reversed1[len1++] = c1;
				reversed2[len2++] = c2;
				if (s1.identity(s2, c1, c2)) {
					identity++;
					similarity++;
				} else if (s1.similarity(s2, c1, c2) > 0) {
					similarity++;
				} else {
				}
				break;
			case Directions.LEFT:
				for (int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = null;
					reversed2[len2++] = --j;
					gaps++;
				}
				break;
			case Directions.STOP:
				stillGoing = false;
			}
			} catch(Exception e) {
				System.out.println(n);
				System.out.println(k+j);
				throw e;
			}
		}

		alignment.setSequence1(reverse(reversed1, len1));
		alignment.setStart1(i);
		alignment.setSequence2(reverse(reversed2, len2));
		alignment.setStart2(j);
		alignment.setIdentity(identity);
		alignment.setGaps(gaps);
		alignment.setSimilarity(similarity);

		return alignment;
	}

	/**
	 * Returns the maximum of 4 double numbers.
	 * 
	 * @param a
	 *            double #1
	 * @param b
	 *            double #2
	 * @param c
	 *            double #3
	 * @param d
	 *            double #4
	 * @return The maximum of a, b, c and d.
	 */
	private static double maximum(double a, double b, double c, double d) {
		if (a > b) {
			if (a > c) {
				return a > d ? a : d;
			} else {
				return c > d ? c : d;
			}
		} else if (b > c) {
			return b > d ? b : d;
		} else {
			return c > d ? c : d;
		}
	}

	/**
	 * Reverses an array of chars
	 * 
	 * @param a
	 * @param len
	 * @return the input array of char reserved
	 */
	private static Integer[] reverse(Integer[] a, int len) {
		Integer[] b = new Integer[len];
		for (int i = len - 1, j = 0; i >= 0; i--, j++) {
			b[j] = a[i];
		}
		return b;
	}
}