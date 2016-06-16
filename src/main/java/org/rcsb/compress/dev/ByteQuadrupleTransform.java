package org.rcsb.compress.dev;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.rcsb.compress.IntegerTransform;

public class ByteQuadrupleTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static int DISTANCE_MIN = 3672;
	private static int DISTANCE_MAX = 3928;
	private static final ByteQuadruples bq = new ByteQuadruples(DISTANCE_MIN, DISTANCE_MAX);
	
//	private static int DISTANCE_MIN = 1800;
//	private static int DISTANCE_MAX = 2000;
	//	private static int DISTANCE_MIN = 3545;
	//	private static int DISTANCE_MAX = 4055;

	public static void main(String[] args) {
		ByteQuadrupleTransform transform = new ByteQuadrupleTransform();
		int[] coords = {2939, 1239, 1900};
		System.out.println("in:  " + Arrays.toString(coords));
		//		int[] results = transform.applyTransformation(coords);
		int[] results = transform.forward(coords);
		System.out.println("forward: " + Arrays.toString(results));
//		int[] rev = transform.reverse(results);
//		System.out.println("reverse: " + Arrays.toString(rev));
//		System.out.println();
	}

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}

	@Override
	public int[] forward(int[] data) {
		int len = data.length / 3;
		int[] out = new int[len*4];

		for (int i = 0; i < len; i++) {
			int[] coords = {data[i],data[i+len], data[i+len+len]};
			int[] results = bq.applyTransformation(coords);

			out[i] = results[0];
			System.out.println(results[0]);
//			if (i > 0) {
//				out[i] = results[0] - out[i-1];
//			}
			out[i + len] = results[1];
			out[i + len + len] = results[2];
			out[i + len + len + len] = results[3];	
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
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
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
