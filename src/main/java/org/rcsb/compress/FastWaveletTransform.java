package org.rcsb.compress;

import java.io.Serializable;
import java.util.Arrays;

import jwave.Transform;
import jwave.TransformBuilder;

public class FastWaveletTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private Transform transform;
	private String waveletName;

	/**
	 * https://en.wikipedia.org/wiki/Discrete_wavelet_transform
	 * adopted to IntegerTransform)
	 */

	public static void main(String[] args) {
		int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532, 868, 586};
	//	int[] data = {11225, 7789, 6238, 2515, 1510, -2120, -389, -209, 435, 804, 708, 1724, 10215, 8596, 5209, 4532};
		System.out.println(Arrays.toString(data));
		
		IntegerTransform t = new AncientEgyptianDecomposition(new FastWaveletTransform("Haar orthogonal"));
		System.out.println(t);
		int[] out = t.forward(data);
		System.out.println(Arrays.toString(out));
		out = t.reverse(out);
		System.out.println(Arrays.toString(out));
	}
	
	
	public FastWaveletTransform(String waveletName) {
		this.transform = TransformBuilder.create("Fast Wavelet Transform",  waveletName);
		this.waveletName = waveletName;
	}
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName() + "(" + waveletName + ")";
	}
	
	@Override
	public int[] forward(int[] input) {
		double[] data = new double[input.length];
		for (int i = 0; i < input.length; i++) {
			data[i] = input[i];
		}
	    
		double[] d = transform.forward(data);
		
		int[] out = new int[input.length];
		for (int i = 0; i < input.length; i++) {
			out[i] = (int) Math.round(d[i]);
		}
//		System.out.println(Arrays.toString(d));
//		System.out.println("Input : " + Arrays.toString(input));
//		System.out.println("Output: " + Arrays.toString(out));
		return out;
	}

	@Override
	public int[] reverse(int[] input) {
		double[] data = new double[input.length];
		for (int i = 0; i < input.length; i++) {
			data[i] = input[i];
		}
	    
		double[] d = transform.reverse(data);
		
		int[] out = new int[input.length];
		for (int i = 0; i < input.length; i++) {
			out[i] = (int) Math.round(d[i]);
		}
		
		return out;
	}
}
