package org.rcsb.compress;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IvqTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static int MAX_ITERATIONS = 20; // 20

	@Override
	public int[] forward(int[] data) {
		int len = data.length/3;

		int[] x = new int[len-1];
		int[] y = new int[len-1];
		int[] z = new int[len-1];
		
		for (int i = 1; i < len+1; i++) {
			x[i-1] = data[i] - data[i-1];
			y[i-1] = data[i+len] - data[i+len-1];
			z[i-1] = data[i+len+len] - data[i+len+len-1];
		}
		
		int[] length = new int[MAX_ITERATIONS];
		int[] v = new int[x.length];
		
		System.out.println("before v: " + Arrays.toString(v));
		System.out.println("before x: " + Arrays.toString(x));
        System.out.println("before y: " + Arrays.toString(y));
        System.out.println("before z: " + Arrays.toString(z));
        
		int iter = 0;
		System.out.print("iter: " + iter);
		int vLength = compress(x, y, z, v);
		System.out.println("after  v: " + Arrays.toString(v));
		System.out.println("after  x: " + Arrays.toString(x));
        System.out.println("after  y: " + Arrays.toString(y));
        System.out.println("after  z: " + Arrays.toString(z));
		length[0] = vLength;
		List<int[]> outs = new ArrayList<>();
		
		while (vLength > 1) {
			vLength = compress(x, y, z, v);
			outs.add(v);
			iter++;
			length[iter]= vLength;
			if (iter == MAX_ITERATIONS-1) break;
		}
		
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
		double[] len = new double[x.length];
		
		int[] xOld = Arrays.copyOf(x, x.length);
		int[] yOld = Arrays.copyOf(y, y.length);
		int[] zOld = Arrays.copyOf(z, z.length);
		
		for (int i = 0; i < x.length; i++) {
			len[i] = Math.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]) + 1;
			
		}
		int length = (int)Math.round(median(len));
		length = (int) (length *1.3);
		
		Arrays.fill(v, 0);
//		double sumSq = 0;

		for (int i = 0; i < x.length; i++) {
				v[i] = CompressUnitVector32Bit.compress3b(x[i]/len[i], y[i]/len[i], z[i]/len[i]);

				double[] d = CompressUnitVector32Bit.decompress3b(v[i]);

//				System.out.println("xold: " + x[i]);
//				int deltaX = (int) Math.round(d[0] * length);
				x[i] -= (int) Math.round(d[0] * length);
//				System.out.println("xcorr: " + (x[i] + deltaX));
				y[i] -= (int) Math.round(d[1] * length);
				z[i] -= (int) Math.round(d[2] * length);

//			}
		} 
		
//		reconsitute(x, y, z, v, length, xOld, yOld, zOld);

		System.out.println("median: " + length);


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
	
	private static double max(double[] len) {
		double max = Double.MIN_VALUE;
		for (double d: len) {
			max = Math.max(d, max);
		}
		return max;
	}
}
