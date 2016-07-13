package org.rcsb.projectm;

import java.io.Serializable;

import javax.vecmath.Point3d;

import org.apache.commons.lang.ArrayUtils;
import org.rcsb.project3.SequenceFeatureInterface;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for the calculation for UsrMomentsFingerprint
 * 
 * @author Chris Li, Michael Wang
 */
public class UsrFeature implements SequenceFeatureInterface<Double>, Serializable {

	private static final long serialVersionUID = 1L;
	private double[] moments;
	private Point3d[] coords;
	// Some setting for calculate similarity score 
	private double match = 1;
	private double mismatch = -1;
	private double gap = -1;
	    
    /**
     * Constructor that will store a double array as a UsrFeature
     * @param moments
     */
	public UsrFeature(double[] moments) {
		this.moments = moments;
	}
	
	/**
	 * Constructor that will store a double array as a UsrFeature and a Point3d array as the 
	 * coordinates of the protein
	 * @param moments
	 * @param coords
	 */
	public UsrFeature(double[] moments, Point3d[] coords) {
		this.moments = moments;
		this.coords = coords;
	}
	
	/**
	 * Constructor that will store a double array as a UsrFeature and update the settings
	 * @param moments
	 * @param gap
	 * @param match
	 * @param mismatch
	 */
	public UsrFeature(double[] moments, double gap, double match, double mismatch) {
		this.moments = moments;
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
		else if (this.get(i).equals(sequence2.get(j)))
			return match;
		else 
			return mismatch;
	}

	@Override
	public boolean identity(SequenceFeatureInterface<Double> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return false;
		}
		// check identity
		else 
			return this.get(i).equals(sequence2.get(j));
	}

	@Override
	public Double get(int index) {
		if (Double.isNaN((double)moments[index]))
			return null;
		else
			return moments[index];
	}

	@Override
	public int length() {
		return moments.length;
	}

	@Override
	public Double[] getSequence() {
		return ArrayUtils.toObject(moments);
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
		return (double)this.get(index);
	}

	@Override
	public Point3d[] getCoords() {
		return coords;
	}
}
