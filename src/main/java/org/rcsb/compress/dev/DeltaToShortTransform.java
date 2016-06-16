package org.rcsb.compress.dev;

import java.io.Serializable;
import java.util.Arrays;

import org.rcsb.compress.IntegerToByteTransform;
import org.rcsb.compress.ShortToByteTransform;

public class DeltaToShortTransform implements IntegerToByteTransform, Serializable {
	private static final long serialVersionUID = 1L;

	
	public static void main(String[] args) {
		int[] data = {0, 100, 1000, 10000, 32000, -32000, -10000, -1000, -100, -0};
		System.out.println(Arrays.toString(data));
		IntegerToByteTransform t = new DeltaToShortTransform();
		byte[] out = t.forward(data);
		System.out.println(Arrays.toString(out));
		int[] back = t.reverse(out);
		System.out.println(Arrays.toString(back));
	}
	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public byte[] forward(int[] data) {
		short[] out = new short[data.length];
//		System.arraycopy(data, 0, out, 0, data.length);
		out[0] = (short) data[0];
		for (int i = out.length-1; i > 0; i--) {
			out[i] = (short) (data[i] - data[i-1]);
		}
		ShortToByteTransform t = new ShortToByteTransformer();
		return t.forward(out);
//		return out;
	}
	
	@Override
	public int[] reverse(byte[] data) {
		ShortToByteTransform t = new ShortToByteTransformer();
		
		short[] delta = t.reverse(data);
	    int[] out = new int[delta.length];
         
	    out[0] = delta[0];
		for (int i = 1; i < delta.length; i++)  {
			out[i] = delta[i-1] + delta[i];
		}
		
		return out;
	}
}
