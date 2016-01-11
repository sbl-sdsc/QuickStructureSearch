package org.rcsb.compress;

import java.io.Serializable;
import java.util.Arrays;

public class CylindricalCoordinateTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static final double scale = 2048;

	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;
		float[] r = new float[len];
		float[] theta = new float[len];
		float[] phi = new float[len];
		int prev = 0;
		for (int i = 0; i < len; i++) { 
			int x = out[i];
			int y = out[i+len];
			r[i] = (float) Math.sqrt(x*x + y*y);
			theta[i] = (float) Math.atan2(y, x);
			out[i] = (int) Math.round(r[i]);
			out[i+len] = (int) Math.round(theta[i]*scale);
			System.out.println(r[i] + " " + out[i] + " " + out[i+len] + " " + out[i+len+len]);
		}
		
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int r = 0;
		int len = data.length/3;
		for (int i = 0; i < len; i++) { 
//			r += out[i];
			r = out[i];
			int theta = out[i+len];
			int phi = out[i+len+len];
			double x = r * Math.sin(theta/scale) * Math.cos(phi/scale);
			double y = r * Math.sin(theta/scale) * Math.sin(phi/scale);
			double z = r * Math.cos(theta/scale);
			out[i] = (int) Math.round(x);
			out[i+len] = (int) Math.round(y);
			out[i+len+len] = (int) Math.round(z);
		}
		
		return out;
	}
}
