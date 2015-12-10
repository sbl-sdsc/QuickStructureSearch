package org.rcsb.compress;

import java.io.Serializable;

public class DeltaTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = out.length-1; i > 0; i--) {
			out[i] = out[i] - out[i-1];
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
	
	private static void printPattern(int[] data) {
		int len = data.length/3;
		if (3*len == data.length) {
		int c = 2000;
		for (int i = 0; i < len; i++) {
			System.out.println(data[i]/c + "\t" + data[i+len]/c + "\t" + data[i+len+len]/c);
		}
		}
	}
	
}
