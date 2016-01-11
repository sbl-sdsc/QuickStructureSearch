package org.rcsb.compress;

import java.io.Serializable;

public class ColToRowTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		int len = data.length/3;
		int n = 0;
		for (int i = 0; i < len; i++) {
			out[n++] = data[i];
			out[n++] = data[i+len];
			out[n++] = data[i+len+len];
		}
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int[] out = new int[data.length];
		
		return out;
	}
}
