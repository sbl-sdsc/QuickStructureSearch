package org.rcsb.hadoop.io;

import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Writable;
/**
 * This class represents a simple biopolymer chain. It supports biopolymers of type: protein, DNA, and RNA chains 
 * and contains the polymer one-letter sequence and the C-alpha atom coordinates of proteins or the
 * P-atom coordinates of nucleic acids.
 * 
 * @author Peter Rose
 *
 */
public class SimplePolymerChain implements Serializable {
	private static final long serialVersionUID = 1L;
	private Writable[] encodedPolymerChain;
	private int start;
	private int end;
	
	public SimplePolymerChain(Writable[] encodedPolymerChain) {
		this.encodedPolymerChain = encodedPolymerChain;
		this.start = SimplePolymerChainCodec.getStartPosition(encodedPolymerChain);
		this.end = SimplePolymerChainCodec.getEndPosition(encodedPolymerChain);
//		System.out.println("New chain: " + start + " - " + end);
	}
	
	/**
	 * @return the polymerType
	 */
	public SimplePolymerType getPolymerType() {
		return SimplePolymerChainCodec.decodePolymerType(encodedPolymerChain);
	}
	/**
	 * @return the coordinates
	 */
	public Point3d[] getCoordinates() {
		return SimplePolymerChainCodec.decodeCoordinates(encodedPolymerChain, start, end);
	}
	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return SimplePolymerChainCodec.decodeSequence(encodedPolymerChain, start, end);
	}
	
	public boolean isProtein() {
		return getPolymerType().equals(SimplePolymerType.PROTEIN);
	}
	
	public boolean isDNA() {
		return getPolymerType().equals(SimplePolymerType.DNA);
	}
	
	public boolean isRNA() {
		return getPolymerType().equals(SimplePolymerType.RNA);
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		sb.append(getPolymerType());
		sb.append(": ");
		sb.append(getSequence());
		sb.append(", ");
		sb.append(Arrays.toString(getCoordinates()));

		return sb.toString();
	}
}
