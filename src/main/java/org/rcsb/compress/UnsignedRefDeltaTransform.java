package org.rcsb.compress;

import java.io.Serializable;

public class UnsignedRefDeltaTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length+1];
		System.arraycopy(data, 0, out, 1, data.length);
		int min = Integer.MAX_VALUE;
		for (int i = out.length-1; i > 1; i--) {
			out[i] = toUnsignedInt(out[i] - out[i-1]);
			min = Math.min(min,  out[i]);
		}
		out[0] = min;
		for (int i = out.length-1; i > 1; i--) {
			out[i] = out[i] - min;
		}
		
		
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int min = data[0];
		int[] out = new int[data.length-1];
		System.arraycopy(data, 1, out, 0, data.length-1);
		for (int i = 1; i < out.length; i++)  {
			out[i] = out[i-1] + toSignedInt(out[i]+min);
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
