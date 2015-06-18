package org.rcsb.project3;

/**
 * This is an interface used to store the fingerprint which is a sequence of feature of the protein
 * @author Chris Li
 *
 * @param <T>
 */
public interface SequenceFeatureInterface<T> {
	/**
	 * Calculate the similarity of a pair of features of the current sequence and the input sequence
	 * @param sequence2	The feature sequence that will be used to compare with
	 * @param i	The feature index of cur sequence
	 * @param j	The feature index of sequence2
	 * @return score for how similar
	 */
	public double similarity(SequenceFeatureInterface<T> sequence2, int i, int j);
	/**
	 * Check the identity of a pair of features of the current sequence and the input sequence
	 * @param sequence2	The feature sequence that will be used to compare with
	 * @param i The feature index of cur sequence
	 * @param j	The feature index of sequence2
	 * @return true if identical, otherwise false
	 */
	public boolean identity(SequenceFeatureInterface<T> sequence2, int i, int j);
	/**
	 * Get the sequence stored in
	 * @return
	 */
	public T[] getSequence();
	/**
	 * Get the feature on the input index of the sequence
	 * @param index
	 * @return
	 */
	public T get(int index);
	/**
	 * Length of the feature sequence
	 * @return
	 */
	public int length();
	/**
	 * Cast a feature at the index to String
	 * @param index
	 * @return
	 */
	public String toString(int index);
}
