package org.rcsb.compress.dev;

import java.io.Serializable;

import org.rcsb.compress.IntegerTransform;

public class SixSphereCoordinateTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static final double scale = 1000000000;

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public int[] forward(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;

		for (int i = 0; i < len; i++) { 
			int x = out[i];
			int y = out[i+len];
			int z = out[i+len+len];
			double dSq = x*x + y*y + z*z;
			double u = x/dSq;
			double v = y/dSq;
			double w = z/dSq;
//			System.out.println(u + " " + v + " " + w);
			out[i] = (int) Math.round(u*scale);
			out[i+len] = (int) Math.round(v*scale);
			out[i+len+len] = (int) Math.round(w*scale);
//			System.out.println(out[i] + " " + out[i+len] + " " + out[i+len+len]);
		}
		
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int[] out = new int[data.length];
		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;

		for (int i = 0; i < len; i++) { 
			double u = out[i]/scale;
			double v = out[i+len]/scale;
			double w = out[i+len+len]/scale;
			double dSq = u*u+ v*v + w*w;
			double x = u/dSq;
			double y = v/dSq;
			double z = w/dSq;
//			System.out.println(u + " " + v + " " + w);
			out[i] = (int) Math.round(x);
			out[i+len] = (int) Math.round(y);
			out[i+len+len] = (int) Math.round(z);
//			System.out.println(out[i] + " " + out[i+len] + " " + out[i+len+len]);
		}
		
		return out;
	}
}
