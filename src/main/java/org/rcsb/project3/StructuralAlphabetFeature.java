package org.rcsb.project3;

import java.io.Serializable;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for the calculation for AngleSequenceFingerprint
 * 
 * @author Chris Li
 */
public class StructuralAlphabetFeature implements SequenceFeatureInterface<String>, Serializable {

	private static final long serialVersionUID = 1L;
	private String[] AlphabetSequence;
	// Some setting for calculate similarity score 
	private double match = 1;
    private double mismatch = -1;
    private double gap = -1;
	
    /**
     * Constructor that will store a double array of angle
     * @param AngleSequence
     */
	public StructuralAlphabetFeature(String[] AlphabetSequence) {
		this.AlphabetSequence = AlphabetSequence;
	}
	
	/**
	 * Constructor that will store a double array of angle and update the settings
	 * @param AngleSequence
	 * @param diff
	 * @param gap
	 * @param match
	 * @param mismatch
	 */
	public StructuralAlphabetFeature(String[] AlphabetSequence, double gap, double match, double mismatch) {
		this.AlphabetSequence = AlphabetSequence;
		this.gap = gap;
		this.match = match;
		this.mismatch = mismatch;
	}
		
	@Override
	public double similarity(SequenceFeatureInterface<String> sequence2, int i, int j) {
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
	public boolean identity(SequenceFeatureInterface<String> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return false;
		}
		// check identity
		else if (this.get(i).equals(sequence2.get(j)))
			return true;
		else 
			return false;
	}

	@Override
	public String get(int index) {
		return AlphabetSequence[index];
	}

	@Override
	public int length() {
		return AlphabetSequence.length;
	}

	@Override
	public String[] getSequence() {
		return AlphabetSequence;
	}

	@Override
	public String toString(int index) {
		return this.get(index);
	}

	@Override
	public double todouble(int index) {
		return 0;
	}
}
