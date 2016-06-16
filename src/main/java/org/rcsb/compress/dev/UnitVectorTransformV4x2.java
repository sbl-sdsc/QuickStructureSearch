package org.rcsb.compress.dev;

import java.io.Serializable;

import org.rcsb.compress.CompressUnitVector8Bit;
import org.rcsb.compress.IntegerTransform;

public class UnitVectorTransformV4x2 implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static int LENGTH = 3800;

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public int[] forward(int[] data) {

//		System.arraycopy(data, 0, out, 0, data.length);
		int len = data.length/3;
//		int[] out = new int[2*len];
		int[] out = new int[5*len];
//		int[] out = new int[len];
//		int deviations = 0;
		
		// calculate average length
		int sum = 0;
		for (int i = 1; i < len; i++) {
			int dx = data[i] - data[i-1];
			int dy = data[i+len] - data[i+len-1];
			int dz = data[i+len+len] - data[i+len+len-1];
			int dSq = dx*dx +dy*dy + dz*dz;
			sum += Math.sqrt(dSq);
		}
		LENGTH = sum/(len-1);
		System.out.println("average length: " + LENGTH);
		
			for (int i = 1; i < len; i++) {
	//			System.out.println(out[i] + "\t" + out[i+len] +"\t" + out[i+len+len]);
				int dx = data[i] - data[i-1];
				int dy = data[i+len] - data[i+len-1];
				int dz = data[i+len+len] - data[i+len+len-1];
//				System.out.println(" Orig: " + dx + "\t" + dy +"\t" + dz);
				int dSq = dx*dx +dy*dy + dz*dz;
				double blen = Math.sqrt(dSq);
				double fx = dx/blen;
				double fy = dy/blen;
				double fz = dz/blen;
//				int deltaD = (int) Math.round(LENGTH - blen);
				int val = CompressUnitVector8Bit.compress3b(fx, fy, fz);
//				System.out.println("orig: " + fx + "\t" + fy +"\t" + fz + "\t" + val);
				double[] d = CompressUnitVector8Bit.decompress3b(val);
//				System.out.println("calc: " + d[0] + "\t" + d[1] +"\t" + d[2] + "\t" + deltaD);
//				int ddx = (int) Math.round(d[0] * (LENGTH - deltaD));
//				int ddy = (int) Math.round(d[1] * (LENGTH - deltaD));
//				int ddz = (int) Math.round(d[2] * (LENGTH - deltaD)); 
				int ddx = (int) Math.round(d[0] * LENGTH);
				int ddy = (int) Math.round(d[1] * LENGTH);
				int ddz = (int) Math.round(d[2] * LENGTH);
//				System.out.println(" Calc: " + ddx + "\t" + ddy +"\t" + ddz);
//				deviations += Math.abs(dx-ddx);
//				deviations += Math.abs(dy-ddy);
//				deviations += Math.abs(dz-ddz);
				System.out.println("delta: " + (dx-ddx) + "\t" + (dy-ddy) +"\t" + (dz-ddz) + "\t val: " + val);
				out[i] = val;
				out[i+len] = dx - ddx;
				out[i+len+len] = dy - ddy;
				out[i+len+len+len] = dz - ddz;
			}
//		System.out.println("deviations: " + deviations);
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
}
