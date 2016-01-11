package org.rcsb.compress;

import java.io.Serializable;
import java.util.Arrays;

import org.apache.commons.lang3.NotImplementedException;

public class PythagoreanTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;
		double[] lengths = new double[len];
		for (int i = 0; i < len; i++) {  //    x^2                y^2                        z^2
			lengths[i] = (int)Math.round(Math.sqrt(out[i]*out[i] + out[i+len]*out[i+len] + out[i+len+len]*out[i+len+len]));
		}
		int median = (int) Math.round(median(lengths));
		System.out.println("median: " + median);
		
		for (int i = 0; i < len; i++) {  //       d^2                y^2                        z^2
			int xe = (int)Math.round(Math.sqrt(median* median - out[i+len]*out[i+len] - out[i+len+len]*out[i+len+len]));
			int dx = Math.abs(out[i]) - xe;
			int udx = toUnsignedInt(dx);
			int sdx = (int)Math.signum(out[i])*udx;
//			System.out.println("x " + out[i] +  " xe: " + xe + " dx: " + dx + " udx: " + udx + " sdx: " + sdx);
			out[i] = sdx;
		}
		
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		throw new NotImplementedException("reverse method not implemented, yet");
	}
	
	private static int toUnsignedInt(final int n) {
	    return (n << 1) ^ (n >> 31);
	}
	
	private static int toSignedInt(final int n) {
	    return (n >>> 1) ^ -(n & 1);
    }
	
	private static double median(double[] len) {
		double[] out = new double[len.length];
		System.arraycopy(len, 0, out, 0, len.length);
		Arrays.sort(out);
		int half = out.length/2;
		if (out.length % 2 == 0) {
			return (out[half-1] + out[half])/2;
		}
		return out[half];
	}
}
