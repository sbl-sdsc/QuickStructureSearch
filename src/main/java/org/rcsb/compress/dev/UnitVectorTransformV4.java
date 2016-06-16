package org.rcsb.compress.dev;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.rcsb.compress.CompressUnitVector8Bit;
import org.rcsb.compress.IntegerTransform;

public class UnitVectorTransformV4 implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public int[] forward(int[] data) {
		int len = data.length/3 - 1;

		int[] x = new int[len];
		int[] y = new int[len];
		int[] z = new int[len];
		
		for (int i = 1; i < len+1; i++) {
			x[i-1] = data[i] - data[i-1];
			y[i-1] = data[i+len] - data[i+len-1];
			z[i-1] = data[i+len+len] - data[i+len+len-1];
		}
		
		int[] o = compress(x, y, z);
		List<int[]> outs = new ArrayList<>();
		
		int iter = 1;
		while (!isZero(o)) {
			System.out.println("iteration: " + iter);
			o = compress(x, y, z);
			iter++;
			if (iter > 20) break;
		}
		
		outs.add(o);
		outs.add(x);
		outs.add(y);
		outs.add(z);
			
		int[] out = new int[outs.size()* o.length];
		for (int i = 0; i < outs.size(); i++) {
			int[] a = outs.get(i);
			System.arraycopy(a, 0, out, i*a.length, a.length);
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
	
	private static int[] compress(int[] x, int[] y, int[] z) {
		double[] len = new double[x.length];

		for (int i = 0; i < x.length; i++) {
			len[i] = (int)Math.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]) + 1;
		}

		double length = (int)median(len);

//		System.out.println("average length: " + length);
		int[] compressedVector = new int[x.length];
		int deviation = 0;
		double sumSq = 0;
		int count = 0;
		for (int i = 0; i < x.length; i++) {
			if (len[i] < 2) continue;
            count +=3;
//			if (Math.max(length, len)/Math.min(length, len) < 1.2) {
				int val = CompressUnitVector8Bit.compress3b(x[i]/len[i], y[i]/len[i], z[i]/len[i]);
				compressedVector[i] = val;

				double[] d = CompressUnitVector8Bit.decompress3b(val);

				x[i] -= (int) Math.round(d[0] * length);
				y[i] -= (int) Math.round(d[1] * length);
				z[i] -= (int) Math.round(d[2] * length);
				//			System.out.println("delta: " + x[i] + ", " + y[i] + ", " + z[i]);
				sumSq = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
				if (Math.abs(x[i]) + Math.abs(y[i]) + Math.abs(z[i]) > 0)
					deviation ++;
			} 
//		}
		System.out.println("rmsd: " + Math.sqrt(sumSq/count));
		System.out.println("dev: " + deviation + " " + (float)deviation/x.length);
		System.out.println("median: " + length);
        System.out.println(Arrays.toString(compressedVector));
        System.out.println(Arrays.toString(x));
        System.out.println(Arrays.toString(y));
        System.out.println(Arrays.toString(z));

		return compressedVector;
	}

	private static boolean isZero(int[] data) {
		for (int d: data) {
			if (d != 0) {
				return false;
			}
		}
		return true;
	}
	
	private static double median(double[] data) {
		double[] out = new double[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		Arrays.sort(out);
		int half = out.length/2;
		if (out.length % 2 == 0) {
			return (out[half-1] + out[half])/2;
		}
		return out[half];
	}
}
