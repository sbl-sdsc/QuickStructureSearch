package org.rcsb.project3;

import java.io.Serializable;

import javax.vecmath.Point3d;

import org.apache.commons.lang.ArrayUtils;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for the calculation for EndToEndDistanceSequenceFingerprint
 * 
 * @author Chris Li
 */
public class EndToEndDistanceSequenceFeature implements SequenceFeatureInterface<Integer>, Serializable {

	private static final long serialVersionUID = 1L;
	private int[] EndToEndSequence;
	private Point3d[] coords;
	// Some setting for calculate similarity score 
	private double match = 1;
	private double mismatch = -1;
	private double gap = -1;
	    
    /**
     * Constructor that will store an int array of DCT1D feature
     * @param DCT1DSequence
     */
	public EndToEndDistanceSequenceFeature(int[] EndToEndSequence) {
		this.EndToEndSequence = EndToEndSequence;
	}
	
	public EndToEndDistanceSequenceFeature(int[] EndToEndSequence, Point3d[] coords) {
		this.EndToEndSequence = EndToEndSequence;
		this.coords = coords;
	}
	
	/**
	 * Constructor that will store an int array of DCT1D feature and update the settings
	 * @param DCT1DSequence
	 * @param gap
	 * @param match
	 * @param mismatch
	 */
	public EndToEndDistanceSequenceFeature(int[] EndToEndSequence, double gap, double match, double mismatch) {
		this.EndToEndSequence = EndToEndSequence;
		this.gap = gap;
		this.match = match;
		this.mismatch = mismatch;
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
		if (Double.isNaN((double)EndToEndSequence[index]))
			return null;
		else
			return EndToEndSequence[index];
	}

	@Override
	public int length() {
		return EndToEndSequence.length;
	}

	@Override
	public Integer[] getSequence() {
		return ArrayUtils.toObject(EndToEndSequence);
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
