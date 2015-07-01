package org.rcsb.hadoop.io;

import java.util.Arrays;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Writable;


public class SimplePolymerChainCodec {
	private static final int SCALE = 1000; 
	private static final double INVERSE_SCALE = 0.001; 
	
	public static Writable[] encodePolymerChain(int type, Point3d[] coords, Integer[] sequence, int gaps) {
		Writable[] writable = new Writable[(coords.length+1)*3 - 2*gaps + 2 + coords.length];
		int n = 0;

		// record polymer type
		writable[n++] = new IntWritable(type);
		
		// record number of points
		writable[n++] = new IntWritable(coords.length);

		int x = 0;
		int y = 0;
		int z = 0;

		// delta encode coordinate values. Gaps in the protein chains
		// are encoded as a maximum integer value.
		for (int i = 0, dx = 0, dy = 0, dz = 0; i < coords.length; i++) {
			Point3d p = coords[i];
			if (p != null) {
				// convert to integers
				x = (int)Math.round(p.x*SCALE);
				y = (int)Math.round(p.y*SCALE);
				z = (int)Math.round(p.z*SCALE);
				
				//delta encode coordinate values as integers
				writable[n++] = new IntWritable(x-dx);
				writable[n++] = new IntWritable(y-dy);
				writable[n++] = new IntWritable(z-dz);
				dx = x;
				dy = y;
				dz = z;
			} else {
				// encode a gap in the protein chain
				writable[n++] = new IntWritable(Integer.MAX_VALUE);
			}
		}

		// record last x,y,z values for validation
		writable[n++] = new IntWritable(x);
		writable[n++] = new IntWritable(y);
		writable[n++] = new IntWritable(z);
		
		// write frame of reference encoded amino acid sequence
		for (Integer s: sequence) {
			writable[n++] = new IntWritable(s);
		}
		
		return writable;
	}

	public static SimplePolymerType decodePolymerType(Writable[] w) {
		int type = ((IntWritable)w[0]).get();
		return SimplePolymerType.values()[type];
	}
	
	public static String decodeSequence(Writable[] encodedPolymerChain, int start, int end) {
		int len = ((IntWritable)encodedPolymerChain[1]).get();	
		
		// skip to position where protein sequence is encoded
		int begin = encodedPolymerChain.length - len;
		
		StringBuilder sb = new StringBuilder();
		for (int i = begin; i < encodedPolymerChain.length; i++) {	
			// amino acid one letter code is encoded as integer that represents its char value
			int charValue = ((IntWritable)encodedPolymerChain[i]).get();	
			sb.append((char)charValue);
		}
		
		String s = sb.toString();
		return s.substring(start, end+1);
	//	return sb.toString();
	}
	
	public static Point3d[] decodeCoordinates(Writable[] encodedPolymerChain, int start, int end) {
		int len = ((IntWritable)encodedPolymerChain[1]).get();
		Point3d[] points = new Point3d[len];
		
		int j = 2;
		int x = 0;
		int y = 0;
		int z = 0;

		for (int i = 0; i < points.length; i++) {
			int v = ((IntWritable)encodedPolymerChain[j++]).get();
			if (v == Integer.MAX_VALUE) {
				points[i] = null; // a gap in the coordinates is represented by a null value
			} else {
				x += v;
				y += ((IntWritable)encodedPolymerChain[j++]).get();
				z += ((IntWritable)encodedPolymerChain[j++]).get();
				points[i] = new Point3d(x*INVERSE_SCALE, y*INVERSE_SCALE, z*INVERSE_SCALE);
			}
		}
		
		// compare the last x, y, z values with the expected values
		if (x != ((IntWritable)encodedPolymerChain[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted");
		}
		if (y != ((IntWritable)encodedPolymerChain[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted");
		}
		if (z != ((IntWritable)encodedPolymerChain[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted");
		}
		
		points = Arrays.copyOfRange(points, start, end+1);
		return points;
	}
	
	/**
	 * Returns the first N-terminal position that has atom coordinates.
	 * 
	 * @param points
	 * @return start position
	 */
	public static int getStartPosition(Writable[] encodedPolymerChain) {
		int len = ((IntWritable)encodedPolymerChain[1]).get();
		
		// find N-terminal start of chain with atom coordinates
		for (int i = 2; i < len+2; i++) {
			int v = ((IntWritable)encodedPolymerChain[i]).get();
			if (v != Integer.MAX_VALUE) {
				return i;
			}
		}
		return 0;
	}
	
	/**
	 * Returns the last (C-terminal) position that has atom coordinates.
	 * 
	 * @param points
	 * @return start position
	 */
	public static int getEndPosition(Writable[] encodedPolymerChain) {
		int len = ((IntWritable)encodedPolymerChain[1]).get();
		int end = 0;

		for (int i = 0, j = 2; i < len; i++) {
			int v = ((IntWritable)encodedPolymerChain[j++]).get();
			if (v != Integer.MAX_VALUE) {
				j+=2;
				end = i;
			}
		}
		
		return end;
	}
}
