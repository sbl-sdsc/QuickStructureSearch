package org.rcsb.project3;

import java.io.Serializable;

import javax.vecmath.Point3d;

import org.apache.commons.lang.ArrayUtils;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for the calculation for AngleSequenceFingerprint
 * 
 * @author Chris Li
 */
public class AngleSequenceFeature implements SequenceFeatureInterface<Double>, Serializable {

	private static final long serialVersionUID = 1L;
	private double[] AngleSequence;
	private Point3d[] coords;
	// Some setting for calculate similarity score 
	private double diff = 3.14/72;
	private double match = 1;
    private double mismatch = -1;
    private double gap = -1;
	
    /**
     * Constructor that will store a double array of angle
     * @param AngleSequence
     */
	public AngleSequenceFeature(double[] AngleSequence) {
		this.AngleSequence = AngleSequence;
	}
	/**
	 * Constructor that also take the protein's 3d coords as an input
	 * @param AngleSequence
	 * @param coords
	 */
	public AngleSequenceFeature(double[] AngleSequence, Point3d[] coords) {
		this.AngleSequence = AngleSequence;
		this.coords = coords;
	}
	
	/**
	 * Constructor that will store a double array of angle and update the settings
	 * @param AngleSequence
	 * @param diff
	 * @param gap
	 * @param match
	 * @param mismatch
	 */
	public AngleSequenceFeature(double[] AngleSequence, double diff, double gap, double match, double mismatch) {
		this.AngleSequence = AngleSequence;
		this.diff = diff;
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
		else if (((this.get(i)+diff) > (double)sequence2.get(j)) && ((this.get(i)-diff) < (double)sequence2.get(j)))
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
		else if (((this.get(i)+diff) > (double)sequence2.get(j)) && ((this.get(i)-diff) < (double)sequence2.get(j)))
			return true;
		else 
			return false;
	}

	@Override
	public Double get(int index) {
		if (Double.isNaN(AngleSequence[index]))
			return null;
		else
			return AngleSequence[index];
	}

	@Override
	public int length() {
		return AngleSequence.length;
	}

	@Override
	public Double[] getSequence() {
		return ArrayUtils.toObject(AngleSequence);
	}

	@Override
	public String toString(int index) {
		if (this.get(index) == null)
			return "NULL";
		else 			
			return String.format("%.2f",this.get(index));
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
