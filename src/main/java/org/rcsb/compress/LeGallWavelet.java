package org.rcsb.compress;

import java.io.Serializable;
import java.util.Arrays;

public class LeGallWavelet implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	/**
	 * Class that represents DiscreteWaveletTransform on Le Gall filter bank
	 * basis using so-called "lifting" scheme. Also it performs
	 * integer-to-integer wavelet transform. It means that for image
	 * data(integer data) decomposition will produce integer coefficients.
	 *
	 * @author Kirilchuk V.E.
	 * @author Peter Rose (adopted to IntegerTransform)
	 */

	public static void main(String[] args) {
		int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532, 868, 586};
	//	int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532};
		System.out.println(Arrays.toString(data));
		
		IntegerTransform t = new AncientEgyptianDecomposition(new DeltaTransform());
		int[] out = t.forward(data);
		out = t.reverse(out);
		System.out.println(Arrays.toString(out));
	}
	public int[] forward(int[] in) {
		int[] out = new int[in.length];
		System.arraycopy(in, 0,  out,  0,  in.length);

		int n = out.length;
		while (n > 1) {
			decompose(out, n);
			n /= 2;
		}
		return out;
	}

	public int[] reverse(int[] in) {
		int[] out = new int[in.length];
		System.arraycopy(in, 0,  out,  0,  in.length);
		
		int n = 1;
		while (n < out.length) {
			reconstruct(out, n);
			n *= 2;
		}
		
		return out;
	}

	private static void reconstruct(int[] x, int n) {
		int endIndex = n * 2;

		// Unpack
		int[] temp = new int[endIndex];
		for (int i = 0; i < n; i++) {
			temp[i * 2] = x[i];
			temp[i * 2 + 1] = x[i + n];
		}
		System.arraycopy(temp, 0, x, 0, endIndex);

		// Undo update
		for (int i = 2; i < endIndex; i += 2) {
			x[i] -= (x[i - 1] + x[i + 1] + 2) >> 2;
		}
		x[0] -= x[1];

		// Undo predict
		for (int i = 1; i < endIndex - 2; i += 2) {
			x[i] += (x[i - 1] + x[i + 1]) >> 1;
		}
		x[endIndex - 1] += x[endIndex - 2];
	}

	/**
	 * +1 level of decomposition.
	 *
	 * <pre>
	 * For example:
	 * We have vector 1,1,1,1,2,2,2,2 where 1 are approximation and 2 are details.
	 * To decompose one more time we need call
	 * decomposeInplace1Lvl([1,1,1,1,2,2,2,2], 4);
	 * 4 - index where details start and approximations ended.
	 * </pre>
	 *
	 * @param x
	 *            vector with approximation and details.
	 * @param endIndex
	 *            index where details start and approximations ended.
	 */
	private static void decompose(int[] x, int endIndex) {
		// Predict
		for (int i = 1; i < endIndex - 2; i += 2) {
			x[i] -= (x[i - 1] + x[i + 1]) >> 1;
		}
		x[endIndex - 1] -= x[endIndex - 2];

		// Update
		for (int i = 2; i < endIndex; i += 2) {
			x[i] += (x[i - 1] + x[i + 1] + 2) >> 2;
		}
		x[0] += x[1];

		// Pack
		int[] temp = new int[endIndex];
		for (int i = 0; i < endIndex; i++) {
			if (i % 2 == 0) {
				temp[i / 2] = x[i];
			} else {
				temp[endIndex / 2 + i / 2] = x[i];
			}
		}

		System.arraycopy(temp, 0, x, 0, endIndex);
	}
}
