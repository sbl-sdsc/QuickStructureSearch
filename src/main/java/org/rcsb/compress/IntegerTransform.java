package org.rcsb.compress;

public interface IntegerTransform {
	public int[] forward(int[] data);
	public int[] reverse(int[] data);
}
