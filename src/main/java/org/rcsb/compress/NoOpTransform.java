package org.rcsb.compress;

import java.io.Serializable;

public class NoOpTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		return forward(data);
	}
	
}
