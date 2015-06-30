package org.rcsb.hadoop.io;

import java.util.Arrays;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Writable;

public class ChainEncoderDecoder {
	private static final int SCALE = 1000; 
	private static boolean truncateTermini = true;
	
	public static Writable[] chainToWritable(int type, Point3d[] coords, Integer[] sequence, int gaps) {
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
	
	public static String writableToSequence(Writable[] w) {
		int len = ((IntWritable)w[1]).get();	
		
		// skip to position where protein sequence is encoded
		int start = w.length - len;
		
		StringBuilder sb = new StringBuilder();
		for (int i = start; i < w.length; i++) {	
			// amino acid one letter code is encoded as integer that represents its char value
			int charValue = ((IntWritable)w[i]).get();	
			sb.append((char)charValue);
		}
		
		return sb.toString();
	}
	
	public static Point3d[] writableToCoordinates(Writable[] w) {
		int len = ((IntWritable)w[1]).get();
		Point3d[] points = new Point3d[len];
		
		int j = 2;
		int x = 0;
		int y = 0;
		int z = 0;

		for (int i = 0; i < points.length; i++) {
			int v = ((IntWritable)w[j++]).get();
			if (v == Integer.MAX_VALUE) {
				points[i] = null; // a gap in the coordinates is represented by a null value
			} else {
				x += v;
				y += ((IntWritable)w[j++]).get();
				z += ((IntWritable)w[j++]).get();
				points[i] = new Point3d(x*SCALE, y*SCALE, z*SCALE);
			}
		}
		
		// compare the last x, y, z values with the expected values
		if (x != ((IntWritable)w[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted");
		}
		if (y != ((IntWritable)w[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted");
		}
		if (z != ((IntWritable)w[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted");
		}
		
		return points;
	}
	
	public static SimplePolymerType writableToPolymerType(Writable[] w) {
		int type = ((IntWritable)w[0]).get();
		return SimplePolymerType.values()[type];
	}
		
	public static SimplePolymerChain writableToSimplePolymerChain(Writable[] w) {
		SimplePolymerType polymerType = writableToPolymerType(w);
		Point3d[] coordinates = writableToCoordinates(w);
		String sequence = writableToSequence(w);
		if (truncateTermini) {
			int start = getStartPosition(coordinates);
			int end = getEndPosition(coordinates);
			coordinates = Arrays.copyOfRange(coordinates, start, end);
			sequence = sequence.substring(start, end+1);
		}
		return new SimplePolymerChain(polymerType, coordinates, sequence);
	}
	
	private static int getStartPosition(Point3d[] points) {
		int start = 0;
		// N-terminal gap (start of chain)
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				start = i;
				break;
			}
		}
		return start;
	}
	
	private static int getEndPosition(Point3d[] points) {
		int end = points.length-1;
		for (int i = points.length-1; i > 0; i--) {
			if (points[i] != null) {
				end = i;
				break;
			}
		}
		return end;
	}
}
