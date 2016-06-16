package org.rcsb.compress.dev;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.rcsb.compress.DeltaTransform;
import org.rcsb.compress.IntegerToByteTransform;
import org.rcsb.compress.IntegerToByteTransformer;
import org.rcsb.compress.IntegerTransform;
import org.rcsb.compress.NoOpTransform;
import org.rcsb.compress.ShortToByteTransform;

public class ByteQuadrupleToByteTransform implements IntegerToByteTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static int DISTANCE_MIN = 3672;
	private static int DISTANCE_MAX = 3928;
	private static final ByteQuadruples bq = new ByteQuadruples(DISTANCE_MIN, DISTANCE_MAX);
	
//	private static int DISTANCE_MIN = 1800;
//	private static int DISTANCE_MAX = 2000;
	//	private static int DISTANCE_MIN = 3545;
	//	private static int DISTANCE_MAX = 4055;

	public static void main(String[] args) {
		ByteQuadrupleToByteTransform transform = new ByteQuadrupleToByteTransform();
		int[] coords = {2939, 1239, 1900};
		System.out.println("in:  " + Arrays.toString(coords));
		//		int[] results = transform.applyTransformation(coords);
//		int[] results = transform.forward(coords);
//		System.out.println("forward: " + Arrays.toString(results));
//		int[] rev = transform.reverse(results);
//		System.out.println("reverse: " + Arrays.toString(rev));
//		System.out.println();
	}

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}

	@Override
	public byte[] forward(int[] in) {
		IntegerTransform t = new DeltaTransform();
		int[] data = t.forward(in);
		int len = data.length / 3;
//		int[] out = new int[len*4];
		short[] indices = new short[len];
		byte[] out = new byte[len*3];

		for (int i = 0; i < len; i++) {
			int[] coords = {data[i],data[i+len], data[i+len+len]};
			int[] results = bq.applyTransformation(coords);

			indices[i] = (short) results[0];
			System.out.println(results[0]);
//			if (i > 0) {
//				out[i] = results[0] - out[i-1];
//			}
			out[i] = (byte) results[1];
			out[i + len] = (byte) results[2];
			out[i + len + len] = (byte) results[3];	
			//			System.out.println("z-corr: " + results[3]);
//			int[] rev = reverse(results);
//			if (rev[0] != coords[0] ||
//					rev[1] != coords[1] ||
//					rev[2] != coords[2] ) {
//				System.out.println("mismatch");
//				System.out.println(Arrays.toString(coords));
//				System.out.println(Arrays.toString(rev));
//			}

		}
		ShortToByteTransform transform = new ShortToByteTransformer();
//		IntegerToByteTransform transform = new IntegerToByteTransformer(new NoOpTransform());
		byte[] bIndices = transform.forward(indices);
		byte[] all = new byte[3*len+ bIndices.length];
		System.arraycopy(out, 0, all, 0, out.length);
		System.arraycopy(bIndices, 0, all, out.length, bIndices.length);
		return all;
	}

	@Override
	public int[] reverse(byte[] data) {
		int len = data.length / 4;
		int[] out = new int[3*len];

		for (int i = 0; i < len; i++) {
			int index = data[i];
			if (index < 0) {
				out[i] = data[i+len];
				out[i + len] = data[i+len+len];
				out[i + len + len] = data[i+len+len+len];
			} else {
				Integer[] quad = null;
				
				// TODO try row order here
				out[i] = quad[1] - data[i+len];
				out[i + len] = quad[2] - data[i + len + len];
				out[i + len + len] = quad[3] - data[i + len+len+len];
			}
		}
		return out;
	}
}
