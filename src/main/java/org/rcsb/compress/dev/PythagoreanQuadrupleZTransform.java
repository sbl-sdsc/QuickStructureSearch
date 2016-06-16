package org.rcsb.compress.dev;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.rcsb.compress.IntegerTransform;

public class PythagoreanQuadrupleZTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private static int DISTANCE_MIN = 3700;
	private static int DISTANCE_MAX = 3900;
	private static int iter = 100;
	private static List<Integer[]> pythagoreanQuadruples;
	static {
//		pythagoreanQuadruples = PythagoreanQuadruples.getReorderedQuadruples(DISTANCE_MIN, DISTANCE_MAX, iter);
		pythagoreanQuadruples = PythagoreanQuadruples.getQuadruples(DISTANCE_MIN, DISTANCE_MAX, iter);
	}
	
	public static void main(String[] args) {
		PythagoreanQuadrupleZTransform transform = new PythagoreanQuadrupleZTransform();
//		Integer[] coords = {3672, 1304, 100};
//		int[] coords = {3600, 1300, 200};
		int[] coords = {2939, 1239, 1900};
		System.out.println("in:  " + Arrays.toString(coords));
//		int[] results = transform.applyTransformation(coords);
		int[] results = transform.forward(coords);
		System.out.println("forward: " + Arrays.toString(results));
		int[] rev = transform.reverse(results);
		System.out.println("reverse: " + Arrays.toString(rev));
		System.out.println();
	}

	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}

	@Override
	public int[] forward(int[] data) {
		int len = data.length / 3;
		int[] out = new int[len*4];
	
		for (int i = 0; i < len; i++) {
			Integer[] coords = 
				{Math.abs(data[i]), 
					Math.abs(data[i + len]), 
					Math.abs(data[i + len + len])};
			int[] results = applyTransformation(coords);

			out[i] = results[0];
			out[i + len] = results[1];
			out[i + len + len] = results[2];
			out[i + len + len + len] = results[3];	
//			System.out.println("z-corr: " + results[3]);
			int[] rev = reverse(results);
			if (rev[0] != coords[0] ||
					rev[1] != coords[1] ||
					rev[2] != coords[2] ) {
				System.out.println("mismatch");
				System.out.println(Arrays.toString(coords));
				System.out.println(Arrays.toString(rev));
			}
					
		}
		return out;
	}

	@Override
	public int[] reverse(int[] data) {
		int len = data.length / 4;
		int[] out = new int[3*len];

		for (int i = 0; i < len; i++) {
			int index = data[i];
			if (index < 0) {
				out[i] = data[i+len];
				out[i + len] = data[i+len+len];
				out[i + len + len] = data[i+len+len+len];
			} else {
			Integer[] quad = pythagoreanQuadruples.get(index);
//			System.out.println("reverse: " + data[i] + "," + data[i+len] + "," + data[i+len+len]);
//			System.out.println("reverse: quad: " + Arrays.toString(quad));
			int x = quad[1] - data[i+len];
			int y = quad[2] - data[i + len + len];
			int dSq = quad[0]*quad[0];
			int xx = x*x;
			int yy = y*y;
			int z = (int)Math.round(Math.sqrt(dSq - xx - yy)) + data[i + len+len+len];
			out[i] = x;
			out[i + len] = y;
			out[i + len + len] = z;
			}
		}
		return out;
	}
	
	private static int[] applyTransformation(Integer[] coords) {
		Arrays.sort(coords, Comparator.reverseOrder());

		int index = -1;
		int dSq = Integer.MAX_VALUE;
		int dxmin = coords[0];
		int dymin = coords[1];
		int dzmin = coords[2];
		
		for (int i = 0; i < pythagoreanQuadruples.size(); i++) {
			Integer[] quad = pythagoreanQuadruples.get(i);
	//		System.out.println(Arrays.toString(quad));
//			if (quad[0] == dist) {
				int dx = quad[1] - coords[0];
				int dy = quad[2] - coords[1];
				int dz = quad[3] - coords[2];
				
				int dTemp = dx*dx + dy*dy + dz*dz;
//				System.out.println(Arrays.toString(quad) + " dSq: " + dTemp);
				if (dTemp < dSq) {
					dSq = dTemp;
					index = i;
					dxmin = dx;
					dymin = dy;
					dzmin = coords[2] - (int) Math.round(Math.sqrt(quad[0]*quad[0] - coords[0]*coords[0] - coords[1]*coords[1]));
				}
			}
//		}
		System.out.println("index: " + index);
		int[] results = {index, dxmin, dymin, dzmin};
		return results;
	}
}
