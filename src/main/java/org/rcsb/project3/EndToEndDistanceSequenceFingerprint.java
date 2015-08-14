package org.rcsb.project3;

import java.io.Serializable;
import java.util.Random;

import javax.vecmath.Point3d;

/**
 * This class generates a fingerprint (signature) for protein chains based upon
 * the end-to-end distance of fragments
 * EndToEndDistance is treated as fragment, and stored as an int array ( for different bin )
 * No similarity for this, may be used for Jaccard Index
 *
 * @author Peter Rose
 */
public class EndToEndDistanceSequenceFingerprint implements SequenceFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	/*
	 * Best default parameter combination: sensitivity: 0.948, specificity: 0.818
	 */
    private int length = 8;
    private double binSize = 3.6943572;
 
    /**
     * Default constructor uses default parameters
     */
    public EndToEndDistanceSequenceFingerprint() {}
    
    /**
     * Constructor with all parameters
     * @param length fragment length
     */
    public EndToEndDistanceSequenceFingerprint(int length, double binSize) {
        this.length = length;
        this.binSize = binSize;
	}
    
    public EndToEndDistanceSequenceFingerprint(int randomSeed) {
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
	public EndToEndDistanceSequenceFeature getFingerprint(Point3d[] coords) {
		int[] features = new int[coords.length-this.length+1];
    	double scale = 1/this.binSize;
    	if (coords.length-this.length < 0) {
    		return new EndToEndDistanceSequenceFeature(features);
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
    		features[i] = (int)Math.round(scale*first.distance(last));
    	}
		return new EndToEndDistanceSequenceFeature(features);
	}
}
