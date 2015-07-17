package org.rcsb.project3;

import java.io.Serializable;

import org.apache.commons.lang.ArrayUtils;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for the calculation for EndToEndDistanceSequenceFingerprint
 * 
 * @author Chris Li
 */
public class EndToEndDistanceDoubleSequenceFeature implements SequenceFeatureInterface<Double>, Serializable {

	private static final long serialVersionUID = 1L;
	private double[] EndToEndSequence;
	// Some setting for calculate similarity score 
	private double match = 1;
	private double mismatch = -1;
	private double gap = -1;
	private double diff = 0.3;
	    
    /**
     * Constructor that will store an int array of DCT1D feature
     * @param DCT1DSequence
     */
	public EndToEndDistanceDoubleSequenceFeature(double[] EndToEndSequence) {
		this.EndToEndSequence = EndToEndSequence;
	}
	
	/**
	 * Constructor that will store an int array of DCT1D feature and update the settings
	 * @param DCT1DSequence
	 * @param gap
	 * @param match
	 * @param mismatch
	 */
	public EndToEndDistanceDoubleSequenceFeature(double[] EndToEndSequence, double gap, double match, double mismatch) {
		this.EndToEndSequence = EndToEndSequence;
		this.gap = gap;
		this.match = match;
		this.mismatch = mismatch;
	}
		
	@Override
	public double similarity(SequenceFeatureInterface<Double> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return gap;
		}
		// check similarity
		else {
			double difference = Math.abs(this.get(i) - sequence2.get(j));
			if (difference < diff)
				return match;
			else {
				double similarity = 1 - difference/1;
				if (similarity > 0)
					return similarity;
				else
					return mismatch;
			}
		}
	}

	@Override
	public boolean identity(SequenceFeatureInterface<Double> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return false;
		}
		// check identity
		else {
			if (this.get(i) + diff > sequence2.get(j) && this.get(i) - diff < sequence2.get(j))
				return true;
			else return false;
		}
	}

	@Override
	public Double get(int index) {
		if (Double.isNaN(EndToEndSequence[index]))
			return null;
		else
			return EndToEndSequence[index];
	}

	@Override
	public int length() {
		return EndToEndSequence.length;
	}

	@Override
	public Double[] getSequence() {
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
		return this.get(index);
	}
}
