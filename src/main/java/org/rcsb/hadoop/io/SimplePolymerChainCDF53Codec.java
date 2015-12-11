package org.rcsb.hadoop.io;

import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.util.Arrays;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Writable;
import org.rcsb.compress.IntegerTransform;

public class SimplePolymerChainCDF53Codec implements Serializable {
	private static final long serialVersionUID = 1L;
	private static final int SCALE = 10; 
	private static final double INVERSE_SCALE = 0.1;
//	private static final int SCALE = 10; 
//	private static final double INVERSE_SCALE = 0.1; 
	
	public static Writable[] encodePolymerChain(int type, Point3d[] caCoords, String sequence, IntegerTransform transform) {	
		int first = getFirstOrderedResidue(caCoords);
		int last = getLastOrderedResidue(caCoords);
		int len = last - first + 1;
		
		int[] ixyz = new int[len*3];
		
		int xOffset = -first;
		int yOffset = xOffset + len;
		int zOffset = yOffset + len;
		
		int x = 0;
		int y = 0;
		int z = 0;
	
		for (int i = first; i < last+1; i++) {
			Point3d p = caCoords[i];
			if (p != null) {
				// convert to integers
				x = (int)Math.round(p.x*SCALE);
				y = (int)Math.round(p.y*SCALE);
				z = (int)Math.round(p.z*SCALE);
				ixyz[i+xOffset] = x;
				ixyz[i+yOffset] = y;
				ixyz[i+zOffset] = z;
			} else {
				ixyz[i+xOffset] = 0; // magic numbers to indicate disordered residue
				ixyz[i+yOffset] = 1;
				ixyz[i+zOffset] = 2;
			}
		}
		
//		int[] txyz = transform.forward(ixyz);
//		int[] rxyz = transform.reverse(txyz);
//		if (ixyz.length != rxyz.length) {
//			System.out.println("Error: orig: " + ixyz.length + " rev: " + rxyz.length);		
//		}
//		for (int i = 0; i < ixyz.length; i++) {
//			if (ixyz[i] != rxyz[i])  System.out.println("Coordinated don't match: " + ixyz[i] + " -> " + rxyz[i]);
//		}
//		System.out.println("Compr: " + ((float)ixyz.length/txyz.length));
		ixyz = transform.forward(ixyz);
		
		int[] seq = new int[sequence.length()];
		for (int i = 0; i < sequence.length(); i++) {
			seq[i] = sequence.charAt(i);
		}
		
//		System.out.println("seq.length: " + seq.length);
		int[] sc = transform.forward(seq);
		int slen = sc.length;
//		System.out.println("sc.length: " + slen);
		
		Writable[] writable = new Writable[6 + ixyz.length + 3 + slen];
		int n = 0;

		// record polymer type
		writable[n++] = new IntWritable(type);
		
		// size of encoded coordinates
		writable[n++] = new IntWritable(ixyz.length);
		
		// record number of points
		writable[n++] = new IntWritable(caCoords.length);
		
		// index of first ordered residue
		writable[n++] = new IntWritable(first);
		
		// index of last ordered residue
		writable[n++] = new IntWritable(last);
		
		// offset to read sequence
		writable[n++] = new IntWritable(6 + ixyz.length + 3);
		
		for (int v: ixyz) {
		   writable[n++] = new IntWritable(v);
		}

		// record last x,y,z values for validation
		writable[n++] = new IntWritable(x);
		writable[n++] = new IntWritable(y);
		writable[n++] = new IntWritable(z);
		
		// write amino acid sequence
		for (int s: sc) {
			writable[n++] = new IntWritable(s);
		}
		
//		String decodedSequence = decodeSequence(writable, transform);
//
//		if (!sequence.equals(decodedSequence)) {
//			System.out.println("Sequence doesn't match: " + sequence.length() + " -> " + decodedSequence.length());
//			for (int i = 0; i < sequence.length(); i++) {
//				if (sequence.charAt(i) != decodedSequence.charAt(i)) System.out.println(sequence.charAt(i) + " -> " + decodedSequence.charAt(i));
//			}
//		    System.out.println(sequence + "$");
//		    System.out.println(decodedSequence + "$");
//		}
		return writable;
	}

	public static SimplePolymerType decodePolymerType(Writable[] w) {
		int type = ((IntWritable)w[0]).get();
		return SimplePolymerType.values()[type];
	}
	
	public static String decodeSequence(Writable[] encodedPolymerChain, IntegerTransform transform) {
		int begin = ((IntWritable)encodedPolymerChain[5]).get();
		int end = encodedPolymerChain.length;
		
		int[] seq = new int[end-begin];
		for (int i = begin; i < encodedPolymerChain.length; i++) {
			seq[i-begin] = ((IntWritable)encodedPolymerChain[i]).get();
		}

		seq = transform.reverse(seq);
		
	    char[] sq = new char[seq.length];
		for (int i = 0; i < sq.length; i++) {
			sq[i] = (char)seq[i];
		}
		
		return new String(sq);
	}
	
	public static Point3d[] decodeCoordinates(Writable[] encodedPolymerChain, IntegerTransform transform) {
		int len = ((IntWritable)encodedPolymerChain[1]).get();
		
		int[] ixyz = new int[len];
	
		int j = 6;
		for (int i = 0; i < ixyz.length; i++) {
			ixyz[i] = ((IntWritable)encodedPolymerChain[j++]).get();
		}

		ixyz = transform.reverse(ixyz);
		
		int cLen = ixyz.length/3;
		int yOffset = cLen;
		int zOffset = yOffset + cLen;
		
		Point3d[] points = new Point3d[cLen];
		
		for (int i = 0; i < cLen; i++) {
			if (ixyz[i] != 0 && ixyz[i+yOffset] != 1 && ixyz[i+zOffset] != 2) {
				points[i] = new Point3d(
						ixyz[i]*INVERSE_SCALE, 
						ixyz[i+yOffset]*INVERSE_SCALE, 
						ixyz[i+zOffset]*INVERSE_SCALE);
			}
		}
		
		// compare the last x, y, z values with the expected values
		if (ixyz[cLen-1] != ((IntWritable)encodedPolymerChain[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted: expected: " + ixyz[cLen-1] + " found " + ((IntWritable)encodedPolymerChain[j-1]).get());
		}
		if (ixyz[cLen-1+yOffset] != ((IntWritable)encodedPolymerChain[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted: expected: " + ixyz[cLen-1+yOffset] + " found " + ((IntWritable)encodedPolymerChain[j-1]).get());
		}
		if (ixyz[cLen-1+zOffset] != ((IntWritable)encodedPolymerChain[j++]).get()) {
			throw new RuntimeException("ERROR: Input file is corrupted: expected: " + ixyz[cLen-1+zOffset] + " found " + ((IntWritable)encodedPolymerChain[j-1]).get());
		}
		
		return points;
	}

	private static int getLastOrderedResidue(Point3d[] points) {
		for (int i = points.length-1; i >=0; i--) {
			if (points[i] != null) {
				return i;
			}
		}
		return points.length;
	}

	private static int getFirstOrderedResidue(Point3d[] points) {
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				return i;
			}
		}
		return points.length;
	}
	
	private static int StringToInt(String string) {
		ByteBuffer buffer = ByteBuffer.wrap(string.getBytes());
		IntBuffer intbuffer = buffer.asIntBuffer();
		intbuffer.rewind();
		
		return buffer.getInt();
	}

}
