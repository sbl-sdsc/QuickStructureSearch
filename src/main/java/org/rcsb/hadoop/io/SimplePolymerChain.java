package org.rcsb.hadoop.io;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Writable;
import org.rcsb.compress.IntegerDeltaZigzagVariableByte;
import org.rcsb.compress.IntegerToByteTransform;
import org.rcsb.compress.dev.DeltaToShortTransform;

public class SimplePolymerChain implements Writable, Serializable {
	private static final long serialVersionUID = 1L;
	private static final int SCALE = 1000; 
	private static final double INVERSE_SCALE = 1.0/SCALE;
	
	private IntegerToByteTransform transform = new IntegerDeltaZigzagVariableByte();
//	private IntegerToByteTransform transform = new DeltaToShortTransform();
	
	private byte polymerType;
	private int coordinateLength;
	private byte[] encodedCoordinates;
	private byte[] encodedSequence;
	
	public SimplePolymerChain() {
	}
	
	public SimplePolymerChain(IntegerToByteTransform transform) {
		this.transform = transform;
	}
	
	public SimplePolymerChain(SimplePolymerChain writable, IntegerToByteTransform transform) {
		this.transform = transform;
		this.polymerType = writable.polymerType;
		this.coordinateLength = writable.coordinateLength;
		this.encodedCoordinates = Arrays.copyOf(writable.encodedCoordinates, writable.encodedCoordinates.length);
		this.encodedSequence = Arrays.copyOf(writable.encodedSequence, writable.encodedSequence.length);
	}
	
	public SimplePolymerChain(SimplePolymerChain writable) {
		this.polymerType = writable.polymerType;
		this.coordinateLength = writable.coordinateLength;
		this.encodedCoordinates = Arrays.copyOf(writable.encodedCoordinates, writable.encodedCoordinates.length);
		this.encodedSequence = Arrays.copyOf(writable.encodedSequence, writable.encodedSequence.length);
	}
	
	public SimplePolymerChain(String sequence, Point3d[] coordinates, int polymerTyper) {
		this.setSequence(sequence);
		this.setCoordinates(coordinates);
		this.setPolymerType(polymerTyper);
	}
	
	public void setPolymerType(int polymerType) {
		this.polymerType = (byte)polymerType;
	}
	
	public int getPolymerType() {
		return this.polymerType;
	}
	
	public boolean isProtein() {
		return getPolymerType() == SimplePolymerType.PROTEIN.ordinal();
	}
	
	public boolean isDNA() {
		return getPolymerType() == SimplePolymerType.DNA.ordinal();
	}
	
	public boolean isRNA() {
		return getPolymerType() == SimplePolymerType.RNA.ordinal();
	}
	
	public void setCoordinates(Point3d[] coordinates) {
		int[] intCoordinates = new int[coordinates.length*3];
		coordinateLength = coordinates.length;
		
		// convert to integers and arrange in x, y, z columns to enable delta encoding
		for (int i = 0; i < coordinates.length; i++) {
			intCoordinates[i] = (int)Math.round(coordinates[i].x*SCALE);
			intCoordinates[i+coordinateLength] = (int)Math.round(coordinates[i].y*SCALE);
			intCoordinates[i+coordinateLength*2] = (int)Math.round(coordinates[i].z*SCALE);
		}
		
//		System.out.println("f: " + Arrays.toString(intCoordinates));
		encodedCoordinates = transform.forward(intCoordinates);
//		System.out.println("e: " + Arrays.toString(encodedCoordinates));
//		int[] outCoordinates = transform.reverse(encodedCoordinates);
//		System.out.println("r: " + Arrays.toString(outCoordinates));
//		
//		if (outCoordinates.length != intCoordinates.length) {
//			throw new Exception("mismatch: " + intCoordinates.length + " - " + outCoordinates.length);
//		}
//		
//		for (int i = 0; i < intCoordinates.length; i++) {
//			if (intCoordinates[i] != outCoordinates[i]) {
//				throw new Exception("coordinate mismatch");
//			}
//		}
	}
	
	public Point3d[] getCoordinates() {
		int[] intCoordinates = transform.reverse(this.encodedCoordinates);
		if (intCoordinates.length != coordinateLength*3) {
			System.err.println("coordinate length inconsistency: " + intCoordinates.length + " - " + coordinateLength);
			System.exit(-1);
		}
		
		Point3d[] coordinates = new Point3d[coordinateLength];
		for (int i = 0; i < coordinateLength; i++) {
		   coordinates[i] = new Point3d(
				   intCoordinates[i]*INVERSE_SCALE,
				   intCoordinates[i+coordinateLength]*INVERSE_SCALE,
				   intCoordinates[i+coordinateLength*2]*INVERSE_SCALE);
		}
		
		return coordinates;		   
	}
	
	public int getCoordinateLength() {
		return coordinateLength;
	}
	
	public void setSequence(String sequence) {
		encodedSequence = sequence.getBytes();
	}
	
	public String getSequence() {
		return new String(encodedSequence);
	}
	
	public int getSequenceLength() {
		return encodedSequence.length;
	}
	
	public int size() {
		return 5 + encodedCoordinates.length + encodedSequence.length;
	}
	
	@Override
	public void write(DataOutput out) throws IOException {
		out.writeByte(this.polymerType);
		out.writeInt(this.coordinateLength);
		out.writeInt(this.encodedCoordinates.length);
		out.write(this.encodedCoordinates);
		out.writeInt(encodedSequence.length);
		out.write(this.encodedSequence);
	}

	@Override
	public void readFields(DataInput in) throws IOException {
		this.polymerType = in.readByte();
		this.coordinateLength = in.readInt();
		
		int encodedCoordinateLength = in.readInt();	
		this.encodedCoordinates = new byte[encodedCoordinateLength];
		in.readFully(this.encodedCoordinates);
		
        int encodedSequenceLength = in.readInt();
        this.encodedSequence = new byte[encodedSequenceLength];
        in.readFully(this.encodedSequence);
	}
}
