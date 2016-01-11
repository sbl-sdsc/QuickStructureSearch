package org.rcsb.compress;

import java.io.Serializable;
import java.nio.ByteBuffer;

/**
 * VariableByte with Delta+Zigzag Encoding
 * Adopted from https://github.com/lemire/JavaFastPFOR
 * 
 * @author MURAOKA Taro http://github.com/koron
 * @author Peter Rose
 */
public class IntegerDeltaZigzagVariableByte implements IntegerToByteTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}

	/**
	 * Transforms and int array into a byte array using Delta + Zigzag + Variable byte encoding.
	 * It applies the following transforms:
	 *    1. Delta encoding (differences between adjacent elements in the input array)
	 *    2. Zigzag encoding to convert deltas to unsigned integers
	 *    3. Variable byte encoding to encode deltas with the minimum number of bytes
	 * 
	 * The first 4 bytes in the byte array encode the length of the original integer array
	 */
	@Override
	public byte[] forward(int[] data) {
		ByteBuffer byteBuf = ByteBuffer.allocate(data.length * 5 + 7);
		byteBuf.putInt(data.length);

		for (int i = 0, p = 0; i < data.length; i++) {

			// delta encoding
			int delta = data[i] - p;
			p = data[i];

			// zigzag encoding (encode as unsigned int)
			delta = (delta << 1) ^ (delta >> 31);

			// variable byte encoding
			switch (Integer.numberOfLeadingZeros(delta)) {
			case 0:
			case 1:
			case 2:
			case 3:
				byteBuf.put((byte) (((delta >>> 28) & 0x7F) | 0x80));
			case 4:
			case 5:
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
				byteBuf.put((byte) (((delta >>> 21) & 0x7F) | 0x80));
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
				byteBuf.put((byte) (((delta >>> 14) & 0x7F) | 0x80));
			case 18:
			case 19:
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
				byteBuf.put((byte) (((delta >>> 7) & 0x7F) | 0x80));
			default:
				byteBuf.put((byte) (delta & 0x7F));
			}
		}

		byteBuf.flip();
		byte[] out = new byte[byteBuf.remaining()];
		byteBuf.get(out);

		return out;
	}

	/**
	 * Transforms a byte array, encoded by Delta + Zigzag + Variable Byte encoding, back into an int array.
	 * It applies the following transforms:
	 *    1. Reverse variable byte encoding
	 *    2. Reverse Zigzag encoding to convert unsigned integers to integers
	 *    3. Reverse delta encoding
	 * 
	 * The first 4 bytes in the byte array encode the length of the original integer array
	 */
	@Override
	public int[] reverse(byte[] data) {
		
		// first 4 bytes are length of int array
		int length = bytes4ToInt(data);
		int[] out = new int[length];

		// decode variable byte and delta+zigzag
		for (int i = 4, j = 0, t = 0, v = 0; i < data.length; i++) { 
			v = (v << 7) + (data[i] & 0x7F);

			if ((data[i] & 0x80) == 0) {
				t = ((v >>> 1) ^ ((v << 31) >> 31)) + t;
				out[j++] = t;
				v = 0;
			}
		}

		return out;
	}
	
	private static int bytes4ToInt(byte[] b) {
		return (b[0] & 0xFF) << 24 | (b[1] & 0xFF) << 16 | (b[2] & 0xFF) << 8 | (b[3]  & 0xFF);
	}
}