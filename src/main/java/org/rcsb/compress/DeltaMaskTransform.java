package org.rcsb.compress;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class DeltaMaskTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public int[] forward(int[] data) {
        int max = 0;
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		for (int i = out.length-1; i > 0; i--) {
			out[i] = out[i] - out[i-1];
			max = Math.max(out[i], max);
		}

		List<Integer> ints = new ArrayList<Integer>();
		int len = out.length/3;
		for (int i = 0; i < len; i++) {
			System.out.println("-----------------------------");
			int[] multipliers = MathToolKit.decompose(Math.abs(out[i])+1);
			System.out.println(Arrays.toString(MathToolKit.decompose(Math.abs(out[i])+1)));
			for (int m: multipliers) {
				ints.add(m);
			}
			multipliers = MathToolKit.decompose(Math.abs(out[i+len])+1);
			System.out.println(Arrays.toString(MathToolKit.decompose(Math.abs(out[i+len])+1)));
			for (int m: multipliers) {
				ints.add(m);
			}
		}
		
		out = new int[ints.size()];
		for (int i = 0; i < out.length; i++) {
			out[i] = ints.get(i);
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
