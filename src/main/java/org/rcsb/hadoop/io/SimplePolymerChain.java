package org.rcsb.hadoop.io;

import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Point3d;
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
	private SimplePolymerType polymerType;
	private Point3d[] coordinates;
	private String sequence;
	
	public SimplePolymerChain(SimplePolymerType polymerType,
			Point3d[] coordinates, String sequence) {
		this.polymerType = polymerType;
		this.coordinates = coordinates;
		this.sequence = sequence;
	}
	
	/**
	 * @return the polymerType
	 */
	public SimplePolymerType getPolymerType() {
		return polymerType;
	}
	/**
	 * @return the coordinates
	 */
	public Point3d[] getCoordinates() {
		return coordinates;
	}
	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}
	
	public boolean isProtein() {
		return polymerType.equals(SimplePolymerType.PROTEIN);
	}
	
	public boolean isDNA() {
		return polymerType.equals(SimplePolymerType.DNA);
	}
	
	public boolean isRNA() {
		return polymerType.equals(SimplePolymerType.RNA);
	}
	
	/**
	 * @param polymerType the polymerType to set
	 */
	public void setPolymerType(SimplePolymerType polymerType) {
		this.polymerType = polymerType;
	}
	/**
	 * @param coordinates the coordinates to set
	 */
	public void setCoordinates(Point3d[] coordinates) {
		this.coordinates = coordinates;
	}
	/**
	 * @param sequence the sequence to set
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(polymerType);
		sb.append(": ");
		sb.append(sequence);
		sb.append(", ");
		sb.append(Arrays.toString(coordinates));

		return sb.toString();
	}
}
