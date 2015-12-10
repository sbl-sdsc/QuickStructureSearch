package org.rcsb.compress;

import java.io.Serializable;

public class Delta4Transform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = out.length-1; i > 1; i--) {
			out[i] = out[i] - out[i-2];
		}
		
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = 4; i < out.length; i++)  {
			out[i] = out[i-4] + out[i];
		}
		
		return out;
	}
	
}
