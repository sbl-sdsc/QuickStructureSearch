package org.rcsb.project3;

import java.io.Serializable;
import java.util.Random;

import javax.vecmath.Point3d;

/**
 * This class generates a fingerprint (signature) for protein chains based upon
 * the end-to-end distance of fragments
 *
 * @author Peter Rose
 */
public class EndToEndDistanceDoubleSequenceFingerprint implements SequenceFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	/*
	 * Best default parameter combination: sensitivity: 0.948, specificity: 0.818
	 */
    private int length = 8;
    private double binSize = 3.6943572;
 
    /**
     * Default constructor uses default parameters
     */
    public EndToEndDistanceDoubleSequenceFingerprint() {}
    
    /**
     * Constructor with all parameters
     * @param length fragment length
     */
    public EndToEndDistanceDoubleSequenceFingerprint(int length, double binSize) {
        this.length = length;
        this.binSize = binSize;
	}
    
    public EndToEndDistanceDoubleSequenceFingerprint(int randomSeed) {
		Random r = new Random(randomSeed);
		this.length = 6 + r.nextInt(3);
		this.binSize = 1.0 + r.nextDouble() * 3;
	}
	
	public float[] getParameters() {
		float[] parameters = {
				this.length,
				(float) this.binSize
				};
		return parameters;
	}

    public String getName() {
    	return this.getClass().getSimpleName() + "_L" + this.length + "B" + this.binSize;
    }

    /**
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
	@Override
	public EndToEndDistanceDoubleSequenceFeature getFingerprint(Point3d[] coords) {
		double[] features = new double[coords.length-this.length+1];
    	double scale = 1/this.binSize;
    	if (coords.length-this.length < 0) {
    		return new EndToEndDistanceDoubleSequenceFeature(features);
    	}

    	for (int i = 0; i < coords.length-this.length+1; i++) {
    		Point3d first = coords[i];
    		Point3d last = coords[i+this.length-1];
    		// skip gaps
    		if (first == null || last  == null) {
    			continue;
    		}
    		// calculate end to end distance of fragment
    		// and bin values
    		features[i] = scale*first.distance(last);
    	}
		return new EndToEndDistanceDoubleSequenceFeature(features);
	}
}
