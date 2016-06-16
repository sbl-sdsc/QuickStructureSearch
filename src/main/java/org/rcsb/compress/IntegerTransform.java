package org.rcsb.compress;

public interface IntegerTransform {
	public String toString();
	public int[] forward(int[] data);
	public int[] reverse(int[] data);
}
