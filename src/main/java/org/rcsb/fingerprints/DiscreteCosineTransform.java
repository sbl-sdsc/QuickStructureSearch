package org.rcsb.fingerprints;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Performs a one dimensional Discrete Cosine Transform on a data set.
 * 
 * @author Alan Yeung, Peter Rose
 */
public class DiscreteCosineTransform implements Serializable {
	private static final long serialVersionUID = 1L;
	private double[][] cosMatrix;
	
	public static void main(String args[]) {
		DiscreteCosineTransform f = new DiscreteCosineTransform();
		double[] data = {10,8,7,6,5,4,3,2,1};
		double[] coefficients = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};
		double[] dct = f.getTransform(data);
		double[] idct = f.getInverseTransform(dct);

		int[] quants = quantize(dct, coefficients);
		double[] vals = dequantize(quants, coefficients, data.length);
		double[] idctq = f.getInverseTransform(vals);
		
		System.out.println("data: " + Arrays.toString(data));
		System.out.println("dct : " + Arrays.toString(dct));
		System.out.println("idct: " + Arrays.toString(idct));
		System.out.println("quant: " + Arrays.toString(quants));
		System.out.println("deqnt: " + Arrays.toString(vals));
		System.out.println("idctq: " + Arrays.toString(idctq));
		System.out.println("rmsd: " + rmsd(data, idctq));
	}
	
	public double[] getTransform(double[] data) {
		if (cosMatrix == null) {
			initializeCosMatrix(data.length);
		}
    	double[] dct = new double[data.length];

    	double sum = 0;
    	for (int m = 0; m < data.length; m++) {
    		sum += data[m];
    	}
    	dct[0] = Math.sqrt(2)/data.length * sum;
    	
    	for (int k = 1; k < data.length; k++){
    		sum = 0;
    		for (int m = 0; m < data.length; m++) {
    			sum += data[m] * cosMatrix[m][k];
    		}
 
    		dct[k] = 2.0/data.length * sum;	
    	}	
    	return dct;
    }
	
	public double[] getInverseTransform(double[] data) {
		if (cosMatrix == null) {
			initializeCosMatrix(data.length);
		}
    	double[] idct = new double[data.length];

    	for (int m = 0; m < data.length; m++) {
    	    double sum = 1/Math.sqrt(2) * data[0];
            for (int k = 1; k < data.length; k++) {
    			sum += data[k] * cosMatrix[m][k];
    		}
            idct[m] = sum;
    	}	
    	return idct;
    }
	
	public static int[] quantize(double[] data, double[] coefficients) {
		int n = Math.min(data.length, coefficients.length);
		int[] values = new int[n];
		for (int i = 0; i < n; i++) {
			values[i] = (int) Math.round(data[i]/coefficients[i]);
			// alternative quantization used in JPEG2000
//			values[i] = (int) (Math.signum(data[i]) * Math.floor(Math.abs(data[i])/coefficients[i]));
		}
		return values;
	}
	
	public static double[] dequantize(int[] quants, double[] coefficients, int length) {
		double[] values = new double[length];
		for (int i = 0; i < coefficients.length; i++) {
			values[i] = quants[i]*coefficients[i];
		}
		return values;
	}
	
	public static double rmsd(double[] x, double[] y) {
		double sum = 0;
		for (int i = 0; i < x.length; i++) {
			sum += (x[i] - y[i]) * (x[i] - y[i]);
		}
		return (Math.sqrt(sum/x.length));
	}
	
	/**
	 * Returns the values for the cos calculations in the DCT summation part.
	 * @return the cos calculations with the x/y value on the horizontal,
	 *    and the u/v values on the vertical
	 */
	private double[][] initializeCosMatrix(int terms) {
		cosMatrix= new double[terms][terms];
		
		for (int m = 0; m < terms; m++) {
			for (int k = 0; k < terms; k++) {
				cosMatrix[m][k] = Math.cos(((2*m + 1) * k*Math.PI)/ (2 * terms));
			}
		}	
		return cosMatrix;	
	}
}
