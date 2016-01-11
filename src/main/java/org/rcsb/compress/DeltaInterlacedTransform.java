package org.rcsb.compress;

import java.io.Serializable;

public class DeltaInterlacedTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = 0; i < data.length-2; i+=2) {
			out[i] = data[i+1];
			out[i+1] = data[i];
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
}
