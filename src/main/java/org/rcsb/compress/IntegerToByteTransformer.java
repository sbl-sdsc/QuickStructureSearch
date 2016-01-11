package org.rcsb.compress;

import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.util.Arrays;

/**
 * 
 * @author Peter Rose
 */
public class IntegerToByteTransformer implements IntegerToByteTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private IntegerTransform transform;

	public static void main(String[] args) {
		int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532, 868, 586};
	//	int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532};
		System.out.println(Arrays.toString(data));
		
		IntegerToByteTransform t = new IntegerToByteTransformer(new SixSphereCoordinateTransform());
		byte[] out = t.forward(data);
		System.out.println(Arrays.toString(out));
		int[] out1 = t.reverse(out);
		System.out.println(Arrays.toString(out1));
	}
	
	public IntegerToByteTransformer(IntegerTransform transform) {
		this.transform = transform;
	}

	/**
	 * Transforms and int array into a byte array
	 */
	@Override
	public byte[] forward(int[] data) {
		int[] data1 = transform.forward(data);
		
		ByteBuffer byteBuf = ByteBuffer.allocate(data1.length * 4);
		for (int value: data1) {
			byteBuf.putInt(value);
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
	public int[] reverse(byte[] data) {
		int[] out = new int[data.length/4];

		ByteBuffer byteBuf = ByteBuffer.wrap(data);
		IntBuffer intBuf = byteBuf.asIntBuffer();
		intBuf.get(out);	
		
		out = transform.reverse(out);
		return out;
	}
}