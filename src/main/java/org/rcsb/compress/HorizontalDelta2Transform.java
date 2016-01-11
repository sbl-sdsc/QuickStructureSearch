package org.rcsb.compress;

import java.io.Serializable;

public class HorizontalDelta2Transform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		int[] temp = new int[3];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;
		for (int i = 1; i < len; i++) {		
			temp[0] = data[i] - data[i-1];
			temp[1] = data[i+len] -  data[i+len-1];
			temp[2] = data[i+len+len] - data[i+len+len-1];
		}
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;
		if (len*3 == data.length) {
			for (int i = 0; i < len; i++) {
				out[i+len] = out[i+len] + out[i];
				out[i+len+len] = out[i+len+len] + out[i];
			}
		}
		
		return out;
	}
}
