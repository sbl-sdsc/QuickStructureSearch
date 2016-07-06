package org.rcsb.projectec;

import java.io.Serializable;

import javax.vecmath.Point3d;

import org.apache.commons.lang.ArrayUtils;
import org.rcsb.project3.SequenceFeatureInterface;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for storing the list of indicies resulting from LibraryFingerprint
 * No similarity for this, may be used for Jaccard Index
 * 
 * @author Emilia Copic
 */
public class LibrarySequenceFeature implements SequenceFeatureInterface<Integer>, Serializable {

	private static final long serialVersionUID = 1L;
	private int[] libSequence;
	private Point3d[] coords;
	// Some setting for calculate similarity score 
	private double match = 1;
	private double mismatch = -1;
	private double gap = -1;
	    
    /**
     * Constructor that will store an int array as EndToEndDistance feature
     * @param EndToEndSequence
     */
	public LibrarySequenceFeature(int[] LibrarySequence) {
		this.libSequence = LibrarySequence;
		
	}
	
	public LibrarySequenceFeature(int[] EndToEndSequence, Point3d[] coords) {
		this.libSequence = EndToEndSequence;
		this.coords = coords;
	}
		
	@Override
	public double similarity(SequenceFeatureInterface<Integer> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return gap;
		}
		// check similarity
		else if (this.get(i).equals(sequence2.get(j)))
			return match;
		else 
			return mismatch;
	}

	@Override
	public boolean identity(SequenceFeatureInterface<Integer> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return false;
		}
		// check identity
		else 
			return this.get(i).equals(sequence2.get(j));
	}

	@Override
	public Integer get(int index) {
		if (Double.isNaN((double)libSequence[index]))
			return null;
		else
			return libSequence[index];
	}

	@Override
	public int length() {
		return libSequence.length;
	}

	@Override
	public Integer[] getSequence() {
		return ArrayUtils.toObject(libSequence);
	}

	@Override
	public String toString(int index) {
		if (this.get(index) == null)
			return "NULL";
		else 			
			return this.get(index).toString();
	}

	@Override
	public double todouble(int index) {
		return (double)(int)this.get(index);
	}

	@Override
	public Point3d[] getCoords() {
		return coords;
	}
}
