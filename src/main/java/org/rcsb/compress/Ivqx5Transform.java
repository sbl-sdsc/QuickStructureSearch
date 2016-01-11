package org.rcsb.compress;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Ivqx5Transform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static int MAX_ITERATIONS = 1; // 20

	@Override
	public int[] forward(int[] data) {
		int len = data.length/3;

		int[] x = new int[len-1];
		int[] y = new int[len-1];
		int[] z = new int[len-1];
		
		int[] xNew = new int[len-1];
		int[] yNew = new int[len-1];
		int[] zNew = new int[len-1];
		
		for (int i = 1; i < len; i++) {
			x[i-1] = data[i] - data[i-1];
			y[i-1] = data[i+len] - data[i+len-1];
			z[i-1] = data[i+len+len] - data[i+len+len-1];
		}
		
		int[] v = new int[x.length];
		
//		System.out.println("before v: " + Arrays.toString(v));
//		System.out.println("before x: " + Arrays.toString(x));
//        System.out.println("before y: " + Arrays.toString(y));
//        System.out.println("before z: " + Arrays.toString(z));
		
		List<int[]> outs = new ArrayList<>();

//			System.out.println("before x: " + Arrays.toString(x));
//	        System.out.println("before y: " + Arrays.toString(y));
//	        System.out.println("before z: " + Arrays.toString(z));
			compress(x, y, z, v);
			outs.add(v);

			System.out.println("v: " + Arrays.toString(v));
			System.out.println("x: " + Arrays.toString(x));
	        System.out.println("y: " + Arrays.toString(y));
	        System.out.println("z: " + Arrays.toString(z));

		
		outs.add(x);
		outs.add(y);
		outs.add(z);
			
		int[] out = new int[outs.size()* v.length];
		int pos = 0;
		for (int i = 0; i < outs.size(); i++) {
			int[] a = outs.get(i);
			System.arraycopy(a, 0, out, pos, a.length);
			pos += a.length;
		}
		
//		System.out.println("iterations: " + outs.size());
//		System.out.println(Arrays.toString(out));


		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;
		if (len*3 == data.length) {
			for (int i = 0; i < len; i++) {
				out[i+len] = out[i+len] + out[i];
				out[i+len+len] = out[i+len+len] + out[i];
			}
		}
		
		return out;
	}
	
	private static int compress(int[] x, int[] y, int[] z, int[] v) {
	    int length = 3800;

	    for (int i = 0; i < x.length; i++) {
	    	v[i] = CompressUnitVector8Bit.compress3b(x[i], y[i], z[i]);

	    	double[] d = CompressUnitVector8Bit.decompress3b(v[i]);

	    	x[i] -= (int) Math.round(d[0] * length);
	    	y[i] -= (int) Math.round(d[1] * length);
	    	z[i] -= (int) Math.round(d[2] * length);
	    } 
		
//		reconsitute(x, y, z, v, length, xOld, yOld, zOld);

//		System.out.println("v: " + Arrays.toString(v));
//		System.out.println("x: " + Arrays.toString(x));
//        System.out.println("y: " + Arrays.toString(y));
//        System.out.println("z: " + Arrays.toString(z));

		return length;
	}

	private static void uncompress(int[] x, int[] y, int[] z, int[] v, int length, int[] xNew, int[] yNew, int[] zNew) {
		for (int i = 0; i < x.length; i++) {
				double[] d = CompressUnitVector8Bit.decompress3b(v[i]);
				xNew[i] = x[i] + (int) Math.round(d[0] * length);
				yNew[i] = y[i] + (int) Math.round(d[1] * length);
				zNew[i] = z[i] + (int) Math.round(d[2] * length);
		} 
	}
	
	private static void reconsitute(int[] x, int[] y, int[] z, int[] v, int length, int[] xOld, int[] yOld, int[] zOld) {
		for (int i = 0; i < x.length; i++) {
				double[] d = CompressUnitVector8Bit.decompress3b(v[i]);
				int xNew = x[i] + (int) Math.round(d[0] * length);
				if (xNew != xOld[i]) 
				System.out.println(xOld[i] + " -> " + xNew);
				int yNew = y[i] + (int) Math.round(d[1] * length);
				if (yNew != yOld[i]) 
				System.out.println(yOld[i] + " -> " + yNew);
				int zNew = z[i] + (int) Math.round(d[2] * length);
				if (zNew != zOld[i]) 
				System.out.println(zOld[i] + " -> " + zNew);
		} 
	}

	private static double median(double[] len) {
		double[] out = new double[len.length];
		System.arraycopy(len, 0, out, 0, len.length);
		Arrays.sort(out);
		int half = out.length/2;
		if (out.length % 2 == 0) {
			return (out[half-1] + out[half])/2;
		}
		return out[half];
	}
	
	private static double percentile(int percentile, double[] len) {
		double[] out = new double[len.length];
		System.arraycopy(len, 0, out, 0, len.length);
		Arrays.sort(out);
		int index = (int) (len.length*percentile*0.01);
		return out[index];
	}
	
	private static double max(double[] len) {
		double max = Double.MIN_VALUE;
		for (double d: len) {
			max = Math.max(d, max);
		}
		return max;
	}
}
