/*
 * $Id: Alignment.java,v 1.11 2006/08/19 00:42:35 ahmed Exp $
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
 * Holds the output of a pairwise sequences alignment.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public final class Alignment<T> {

	/**
	 * Gap character
	 */
	public static final char GAP = '-';

	/**
	 * Gap open cost
	 */
	private double open;

	/**
	 * Gap extend cost
	 */
	private double extend;

	/**
	 * Alignment score
	 */
	private double score;

	/**
	 * Aligned sequence #1
	 */
	private Integer[] sequence1;

	/**
	 * Alignment start location in sequence #1
	 */
	private int start1;

	/**
	 * Aligned sequence #2
	 */
	private Integer[] sequence2;

	/**
	 * Alignment start location in sequence #2
	 */
	private int start2;

	/**
	 * Count of identical locations
	 */
	private int identity;

	/**
	 * Count of similar locations
	 */
	private int similarity;

	/**
	 * Count of gap locations
	 */
	private int gaps;

	private SequenceFeatureInterface<T> originalSequence1;

	private SequenceFeatureInterface<T> originalSequence2;

	/**
	 * Constructor for Alignment
	 */

	public Alignment() {
		super();
	}

	/**
	 * @return Returns the extend.
	 */
	public double getExtend() {
		return extend;
	}

	/**
	 * @param extend
	 *            The extend to set.
	 */
	public void setExtend(double extend) {
		this.extend = extend;
	}

	/**
	 * @return Returns the open.
	 */
	public double getOpen() {
		return open;
	}

	/**
	 * @param open
	 *            The open to set.
	 */
	public void setOpen(double open) {
		this.open = open;
	}

	/**
	 * @return Returns the score.
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score
	 *            The score to set.
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * Returns the length of the alignment
	 * 
	 * @return Alignment length
	 */
	public int getLength() {
		return this.sequence1.length;
	}

	/**
	 * @return Returns the sequence1.
	 */
	public Integer[] getSequence1() {
		return sequence1;
	}

	/**
	 * @param sequence1
	 *            The sequence1 to set.
	 */
	public void setSequence1(Integer[] sequence1) {
		this.sequence1 = sequence1;
	}

	/**
	 * @return Returns the sequence2.
	 */
	public Integer[] getSequence2() {
		return sequence2;
	}

	/**
	 * @param sequence2
	 *            The sequence2 to set.
	 */
	public void setSequence2(Integer[] sequence2) {
		this.sequence2 = sequence2;
	}

	/**
	 * @return Returns the start1.
	 */
	public int getStart1() {
		return start1;
	}

	/**
	 * @param start1
	 *            The start1 to set.
	 */
	public void setStart1(int start1) {
		this.start1 = start1;
	}

	/**
	 * @return Returns the start2.
	 */
	public int getStart2() {
		return start2;
	}

	/**
	 * @param start2
	 *            The start2 to set.
	 */
	public void setStart2(int start2) {
		this.start2 = start2;
	}

	/**
	 * @return Returns the gaps.
	 */
	public int getGaps() {
		return gaps;
	}

	/**
	 * @param gaps
	 *            The gaps to set.
	 */
	public void setGaps(int gaps) {
		this.gaps = gaps;
	}

	/**
	 * @return Returns the identity.
	 */
	public int getIdentity() {
		return identity;
	}

	/**
	 * @param identity
	 *            The identity to set.
	 */
	public void setIdentity(int identity) {
		this.identity = identity;
	}

	/**
	 * @return Returns the similarity.
	 */
	public int getSimilarity() {
		return similarity;
	}

	/**
	 * @param similarity
	 *            The similarity to set.
	 */
	public void setSimilarity(int similarity) {
		this.similarity = similarity;
	}

	/**
	 * Calculate the score of the alignment, not using the score field (the
	 * function only uses sequence1, sequence2, matrix and gap penalties).
	 * 
	 * @return the calculated score (By: Bram Minnaert)
	 */
	public double calculateScore() {
		// The calculated score
		double calcScore = 0;

		// In the previous step there was a gap in the first sequence
		boolean previous1wasGap = false;

		// In the previous step there was a gap in the second sequence
		boolean previous2wasGap = false;

		int start = 0;
		int end = sequence1.length - 1;

		Integer c1, c2; // the next character
		for (int i = start; i <= end; i++) {
			c1 = sequence1[i];
			c2 = sequence2[i];
			// the next character in the first sequence is a gap
			if (c1 == null) {
				if (previous1wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = true;
				previous2wasGap = false;
			}
			// the next character in the second sequence is a gap
			else if (c2 == null) {
				if (previous2wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = false;
				previous2wasGap = true;
			}
			// the next characters in boths sequences are not gaps
			else {
				calcScore += originalSequence1.similarity(originalSequence2, c1, c2);
				previous1wasGap = false;
				previous2wasGap = false;
			}
		}
		return calcScore;
	}

	/**
	 * Calculate the score of the alignment without the terminal gaps.
	 * 
	 */
	public double getScoreWithNoTerminalGaps() {
		// The calculated score
		double calcScore = 0;

		// In the previous step there was a gap in the first sequence
		boolean previous1wasGap = false;

		// In the previous step there was a gap in the second sequence
		boolean previous2wasGap = false;

		int start = 0;
		int end = sequence1.length - 1;

		if (sequence1[start] == GAP) {
			while (sequence1[start] == GAP) {
				start++;
			}
		} else if (sequence2[start] == GAP) {
			while (sequence2[start] == GAP) {
				start++;
			}
		}

		if (sequence1[end] == GAP) {
			while (sequence1[end] == GAP) {
				end--;
			}
		} else if (sequence2[end] == GAP) {
			while (sequence2[end] == GAP) {
				end--;
			}
		}

		Integer c1, c2; // the next character
		for (int i = start; i <= end; i++) {
			c1 = sequence1[i];
			c2 = sequence2[i];
			// the next character in the first sequence is a gap
			if (c1 == GAP) {
				if (previous1wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = true;
				previous2wasGap = false;
			}
			// the next character in the second sequence is a gap
			else if (c2 == GAP) {
				if (previous2wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = false;
				previous2wasGap = true;
			}
			// the next characters in boths sequences are not gaps
			else {
				calcScore += originalSequence1.similarity(originalSequence2, c1, c2);
				previous1wasGap = false;
				previous2wasGap = false;
			}
		}
		return calcScore;
	}

	/**
	 * Check if the calculated score matches the field score.
	 * 
	 * @return true if equal, else false. (By: Bram Minnaert)
	 */
	public boolean checkScore() {
		if (calculateScore() == score) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Returns original {@link Sequence} #1
	 * 
	 * @return original {@link Sequence} #1
	 */
	public SequenceFeatureInterface<T> getOriginalSequence1() {
		return originalSequence1;
	}

	/**
	 * 
	 * @param originalSequence1
	 */
	public void setOriginalSequence1(SequenceFeatureInterface<T> originalSequence1) {
		this.originalSequence1 = originalSequence1;
	}

	/**
	 * Returns original {@link Sequence} #2
	 * 
	 * @return original {@link Sequence} #2
	 */
	public SequenceFeatureInterface<T> getOriginalSequence2() {
		return originalSequence2;
	}

	/**
	 * 
	 * @param originalSequence2
	 */
	public void setOriginalSequence2(SequenceFeatureInterface<T> originalSequence2) {
		this.originalSequence2 = originalSequence2;
	}

	/**
	 * Returns the number of gaps of the aligned sequence #1
	 * 
	 * @return the number of gaps of the aligned sequence #1
	 */
	public int getGaps1() {
		int count = 0;
		for (int i = 0, n = sequence1.length; i < n; i++) {
			if (sequence1[i] == Alignment.GAP) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Returns the number of gaps of the aligned sequence #2
	 * 
	 * @return the number of gaps of the aligned sequence #2
	 */
	public int getGaps2() {
		int count = 0;
		for (int i = 0, n = sequence2.length; i < n; i++) {
			if (sequence2[i] == Alignment.GAP) {
				count++;
			}
		}
		return count;
	}
}