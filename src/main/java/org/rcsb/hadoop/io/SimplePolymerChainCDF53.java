package org.rcsb.hadoop.io;

import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.Writable;
import org.rcsb.compress.IntegerTransform;
/**
 * This class represents a simple biopolymer chain. It supports biopolymers of type: protein, DNA, and RNA chains 
 * and contains the polymer one-letter sequence and the C-alpha atom coordinates of proteins or the
 * P-atom coordinates of nucleic acids.
 * 
 * @author Peter Rose
 *
 */
public class SimplePolymerChainCDF53 implements Serializable {
	private static final long serialVersionUID = 1L;
	private Writable[] encodedPolymerChain;
	private IntegerTransform transform;
	
	public SimplePolymerChainCDF53(Writable[] encodedPolymerChain, IntegerTransform transform) {
		this.encodedPolymerChain = encodedPolymerChain;
		this.transform = transform;
	}
	
	/**
	 * @return the polymerType
	 */
	public SimplePolymerType getPolymerType() {
		return SimplePolymerChainCDF53Codec.decodePolymerType(encodedPolymerChain);
	}
	
	/**
	 * @return size of internal data structures
	 */
	public int size() {
		return encodedPolymerChain.length;
	}
	
	/**
	 * @return the coordinates
	 */
	public Point3d[] getCoordinates() {
		return SimplePolymerChainCDF53Codec.decodeCoordinates(encodedPolymerChain, transform);
	}
	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return SimplePolymerChainCDF53Codec.decodeSequence(encodedPolymerChain, transform);
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
