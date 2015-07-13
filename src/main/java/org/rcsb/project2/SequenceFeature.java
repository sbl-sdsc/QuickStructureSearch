package org.rcsb.project2;


public interface SequenceFeature<T extends SequenceFeature<T, V>, V> {
	public double similarity(SequenceFeature<T, V> sequence2, int i, int j);

	public boolean identity(SequenceFeature<T, V> sequence2, int i, int j);

	public V[] getSequence();

	public V get(int index);

	public int length();

	public String toString(int index);
}
