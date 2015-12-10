package org.rcsb.compress;

import java.io.Serializable;

public class PythagoreanTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;
		for (int i = 0; i < len; i++) {  //       d^2                y^2                        z^2
			int xe = (int)Math.round(Math.sqrt(3800* 3800 - out[i+len]*out[i+len] - out[i+len+len]*out[i+len+len]));
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
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = 1; i < out.length; i++)  {
			out[i] = out[i-1] + out[i];
		}
		
		return out;
	}
	
	private static int toUnsignedInt(final int n) {
	    return (n << 1) ^ (n >> 31);
	}
	
	private static int toSignedInt(final int n) {
	    return (n >>> 1) ^ -(n & 1);
    }
}
