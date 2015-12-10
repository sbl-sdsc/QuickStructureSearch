package org.rcsb.compress;

import java.io.Serializable;
import java.util.Arrays;

public class HarrWavelet implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	/**
	 * https://en.wikipedia.org/wiki/Discrete_wavelet_transform
	 * adopted to IntegerTransform)
	 */

	public static void main(String[] args) {
		int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532, 868, 586};
	//	int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532};
		System.out.println(Arrays.toString(data));
		
		IntegerTransform t = new AncientEgyptianDecomposition(new HarrWavelet());
		int[] out = t.forward(data);
		System.out.println(Arrays.toString(out));
		out = t.reverse(out);
		System.out.println(Arrays.toString(out));
	}
	
	public int[] forward(int[] input) {
	    // This function assumes that input.length=2^n, n>1
	    int[] output = new int[input.length];
	    if (input.length == 1) {
	    	output[0] = input[0];
	    	return output;
	    }

	    for (int length = input.length >> 1; ; length >>= 1) {
	        // length = input.length / 2^n, WITH n INCREASING to log_2(input.length)
	        for (int i = 0; i < length; ++i) {
	            int sum = input[i * 2] + input[i * 2 + 1];
	            int difference = input[i * 2] - input[i * 2 + 1];
	            output[i] = sum;
	            output[length + i] = difference;
	        }
	        if (length == 1) {
	            return output;
	        }

	        //Swap arrays to do next iteration
	        System.arraycopy(output, 0, input, 0, length << 1);
	    }
	}

	public int[] reverse(int[] in) {
		return forward(in);// TODO implement reverse
	}
}
