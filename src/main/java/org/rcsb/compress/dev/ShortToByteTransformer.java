package org.rcsb.compress.dev;

import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.nio.ShortBuffer;
import java.util.Arrays;

import org.rcsb.compress.NoOpTransform;
import org.rcsb.compress.ShortToByteTransform;
import org.rcsb.compress.dev.SixSphereCoordinateTransform;

/**
 * 
 * @author Peter Rose
 */
public class ShortToByteTransformer implements ShortToByteTransform, Serializable {
	private static final long serialVersionUID = 1L;

	public static void main(String[] args) {
		short[] data = {10000, 11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532, 868, 586};
	//	int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532};
		System.out.println(Arrays.toString(data));
		
		ShortToByteTransform t = new ShortToByteTransformer();
		byte[] out = t.forward(data);
		System.out.println(Arrays.toString(out));
		short[] out1 = t.reverse(out);
		System.out.println(Arrays.toString(out1));
	}
	
	public ShortToByteTransformer() {
	}
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}

	/**
	 * Transforms and int array into a byte array
	 */
	@Override
	public byte[] forward(short[] data) {
		
		ByteBuffer byteBuf = ByteBuffer.allocate(data.length * 2);
		for (short value: data) {
			byteBuf.putShort(value);
		}

		byteBuf.flip();
		byte[] out = new byte[byteBuf.remaining()];
		byteBuf.get(out);

		return out;
	}

	/**
	 * Transforms a byte array into an integer array
	 */
	@Override
	public short[] reverse(byte[] data) {
		short[] out = new short[data.length/2];

		ByteBuffer byteBuf = ByteBuffer.wrap(data);
		ShortBuffer intBuf = byteBuf.asShortBuffer();
		intBuf.get(out);	

		return out;
	}
}