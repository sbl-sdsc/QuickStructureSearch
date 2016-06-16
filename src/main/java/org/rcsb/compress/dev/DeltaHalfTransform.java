package org.rcsb.compress.dev;

import java.io.Serializable;
import java.util.Arrays;

import org.rcsb.compress.IntegerTransform;

public class DeltaHalfTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = out.length-1; i > 0; i--) {
			out[i] = (out[i] - out[i-1])/2;
		}
		return out;
	}
	
	@Override
	public int[] reverse(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = 1; i < out.length; i++)  {
			out[i] = (out[i-1] + out[i])*2;
		}
		
		return out;
	}
}
