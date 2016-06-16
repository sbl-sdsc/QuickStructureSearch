package org.rcsb.compress.dev;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.rcsb.compress.CompressUnitVector8Bit;
import org.rcsb.compress.IntegerTransform;

public class Ivqx5Transform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static int MAX_ITERATIONS = 1; // 20

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}

	@Override
	public int[] forward(int[] data) {
		int len = data.length / 3;

		int[] x = new int[len];
		int[] y = new int[len];
		int[] z = new int[len];
		x[0] = data[0];
		y[0] = data[len];
		z[0] = data[len+len];

		for (int i = 1; i < len; i++) {
			x[i] = data[i] - data[i-1];
			y[i] = data[i + len] - data[i + len-1];
			z[i] = data[i + len + len] - data[i + len + len-1];
		}

		int[] v = new int[x.length];

		List<int[]> outs = new ArrayList<>();

		System.out.println("before x: " + Arrays.toString(x));
		System.out.println("before y: " + Arrays.toString(y));
		System.out.println("before z: " + Arrays.toString(z));
		int[] length = new int[1];
		length[0] = compress(x, y, z, v);
		System.out.println("compressed median length: " + length[0]);

		System.out.println("v: " + Arrays.toString(v));
		System.out.println("x: " + Arrays.toString(x));
		System.out.println("y: " + Arrays.toString(y));
		System.out.println("z: " + Arrays.toString(z));

		outs.add(x);
		outs.add(y);
		outs.add(z);
		outs.add(v);
		outs.add(length);

		int[] out = new int[4 * x.length + 1];
		int pos = 0;
		for (int i = 0; i < outs.size(); i++) {
			int[] a = outs.get(i);
			System.arraycopy(a, 0, out, pos, a.length);
			pos += a.length;
		}

		System.out.println("outs.size: " + outs.size());
		System.out.println("outs: " + Arrays.toString(out));

		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int length = (data.length-1)/4;
		int[] x = new int[length];
		int[] y = new int[length];
		int[] z = new int[length];
		int[] v = new int[length];

		x[0] = data[0];
		y[0] = data[length];
		z[0] = data[2* length];
		
		for (int i = 1; i < length; i++) {
			x[i] = data[i];
			y[i] = data[i+length];
			z[i] = data[i+length+length];
			v[i] = data[i + 3*length];
		}

		int vectorLength = data[data.length-1];
		System.out.println("median length: " + vectorLength);
		for (int i = 0; i < x.length; i++) {
			double[] d = CompressUnitVector8Bit.decompress3b(v[i]);

			x[i] += (int) Math.round(d[0] * vectorLength);
			y[i] += (int) Math.round(d[1] * vectorLength);
			z[i] += (int) Math.round(d[2] * vectorLength);
		}
		
		int[] out = new int[3*length];
		for (int i = 0; i < length; i++) {
			out[i] = x[i];
			out[i+length] = y[i];
			out[i+length+length] = z[i];
		}

		return out;
	}

	private static int compress(int[] x, int[] y, int[] z, int[] v) {
		int length = getMedianLength(x, y, z);

		for (int i = 0; i < x.length; i++) {
			v[i] = CompressUnitVector8Bit.compress3b(x[i], y[i], z[i]);

			double[] d = CompressUnitVector8Bit.decompress3b(v[i]);

			x[i] -= (int) Math.round(d[0] * length);
			y[i] -= (int) Math.round(d[1] * length);
			z[i] -= (int) Math.round(d[2] * length);
		}

//		uncompress(x, y, z, v, length);
//
//		System.out.println("v: " + Arrays.toString(v));
//		System.out.println("x: " + Arrays.toString(x));
//		System.out.println("y: " + Arrays.toString(y));
//		System.out.println("z: " + Arrays.toString(z));

		return length;
	}

	private static void uncompress(int[] x, int[] y, int[] z, int[] v,
			int length) {
		for (int i = 0; i < x.length; i++) {
			double[] d = CompressUnitVector8Bit.decompress3b(v[i]);
			x[i] = x[i] + (int) Math.round(d[0] * length);
			y[i] = y[i] + (int) Math.round(d[1] * length);
			z[i] = z[i] + (int) Math.round(d[2] * length);
		}
	}

	private static void reconsitute(int[] x, int[] y, int[] z, int[] v,
			int length, int[] xOld, int[] yOld, int[] zOld) {
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

	private static int getMedianLength(int[] x, int[] y, int[] z) {
		double[] lengths = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			lengths[i] = Math.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
		}
		return (int)Math.round(median(lengths));
	}
	
	private static double median(double[] len) {
		double[] out = new double[len.length];
		System.arraycopy(len, 0, out, 0, len.length);
		Arrays.sort(out);
		int half = out.length / 2;
		if (out.length % 2 == 0) {
			return (out[half - 1] + out[half]) / 2;
		}
		return out[half];
	}

	private static double percentile(int percentile, double[] len) {
		double[] out = new double[len.length];
		System.arraycopy(len, 0, out, 0, len.length);
		Arrays.sort(out);
		int index = (int) (len.length * percentile * 0.01);
		return out[index];
	}

	private static double max(double[] len) {
		double max = Double.MIN_VALUE;
		for (double d : len) {
			max = Math.max(d, max);
		}
		return max;
	}
}
