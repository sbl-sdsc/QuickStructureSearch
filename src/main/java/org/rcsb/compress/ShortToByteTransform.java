package org.rcsb.compress;

/**
 * 
 * @author Peter Rose
 *
 */
public interface ShortToByteTransform {
	public String toString();
	public byte[] forward(short[] data);
	public short[] reverse(byte[] data);
}
