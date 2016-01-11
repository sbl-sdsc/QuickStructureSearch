package org.rcsb.compress;

/**
 * 
 * @author Peter Rose
 *
 */
public interface IntegerToByteTransform {
	public String toString();
	public byte[] forward(int[] data);
	public int[] reverse(byte[] data);
}
