package org.rcsb.project3;

import org.apache.commons.lang.ArrayUtils;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for the calculation for DCT1DSequenceFingerprint
 * 
 * @author Chris Li
 */
public class DCT1DSequenceFeature implements SequenceFeatureInterface<Integer> {

	private int[] DCT1DSequence;
	// Some setting for calculate similarity score 
	private double match = 1;
	private double mismatch = -1;
	private double gap = -1;
	    
    /**
     * Constructor that will store a int array of DCT1D feature
     * @param DCT1DSequence
     */
	public DCT1DSequenceFeature(int[] DCT1DSequence) {
		this.DCT1DSequence = DCT1DSequence;
	}
	
	/**
	 * Constructor that will store a int array of DCT1D feature and update the settings
	 * @param DCT1DSequence
	 * @param gap
	 * @param match
	 * @param mismatch
	 */
	public DCT1DSequenceFeature(int[] DCT1DSequence, double gap, double match, double mismatch) {
		this.DCT1DSequence = DCT1DSequence;
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
		if (Double.isNaN((double)DCT1DSequence[index]))
			return null;
		else
			return DCT1DSequence[index];
	}

	@Override
	public int length() {
		return DCT1DSequence.length;
	}

	@Override
	public Integer[] getSequence() {
		return ArrayUtils.toObject(DCT1DSequence);
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
}
