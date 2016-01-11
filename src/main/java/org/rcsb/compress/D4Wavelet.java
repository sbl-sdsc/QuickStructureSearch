package org.rcsb.compress;

import java.io.Serializable;
import java.util.Arrays;

public class D4Wavelet implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	/**
	 * https://en.wikipedia.org/wiki/Discrete_wavelet_transform
	 * adopted to IntegerTransform)
	 */

	public static void main(String[] args) {
		int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532, 868, 586};
	//	int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532};
		System.out.println(Arrays.toString(data));
		
		IntegerTransform t = new AncientEgyptianDecomposition(new D4Wavelet());
		int[] out = t.forward(data);
		System.out.println(Arrays.toString(out));
		out = t.reverse(out);
		System.out.println(Arrays.toString(out));
	}
	
	public int[] forward(int[] input) {
	    Daub d4 = new Daub();
	    double[] data = new double[input.length];
	    for (int i = 0; i < input.length; i++) {
	    	data[i] = input[i];
	    }
	    	
	    d4.daubTrans(data);
	    
	    int[] out = new int[input.length];
	    for (int i = 0; i < input.length; i++) {
	    	out[i] = (int) Math.round(data[i]);
	    }
	    return out;
	}

	public int[] reverse(int[] input) {
		Daub d4 = new Daub();
	    double[] data = new double[input.length];
	    for (int i = 0; i < input.length; i++) {
	    	data[i] = input[i];
	    }
	    	
	    d4.invDaubTrans(data);
	    
	    int[] out = new int[input.length];
	    for (int i = 0; i < input.length; i++) {
	    	out[i] = (int) Math.round(data[i]);
	    }
	    return out;
	}
}
