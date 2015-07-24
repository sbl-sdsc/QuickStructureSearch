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
public class RmsdDistanceSequenceFingerprint implements SequenceFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	/*
	 * Best default parameter combination: sensitivity: 0.948, specificity: 0.818
	 */
    private int length = 8;
 
    /**
     * Default constructor uses default parameters
     */
    public RmsdDistanceSequenceFingerprint() {}
    
    /**
     * Constructor with all parameters
     * @param length fragment length
     */
    public RmsdDistanceSequenceFingerprint(int length, double binSize) {
        this.length = length;
	}
    
    public RmsdDistanceSequenceFingerprint(int randomSeed) {
		Random r = new Random(randomSeed);
		this.length = 6 + r.nextInt(3);
	}
	
	public float[] getParameters() {
		float[] parameters = {
				this.length
				};
		return parameters;
	}

    public String getName() {
    	return this.getClass().getSimpleName() + "_L" + this.length;
    }

    /**
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
	@Override
	public EndToEndDistanceDoubleSequenceFeature getFingerprint(Point3d[] coords) {
		double[] features = new double[coords.length-this.length+1];
    	if (coords.length-this.length < 0) {
    		return new EndToEndDistanceDoubleSequenceFeature(features);
    	}

    	for (int i = 0; i < coords.length-this.length+1; i++) {
    		features[i] = 0;
    		for (int j = i; j < i + this.length - 1; j++) {
        		Point3d p1 = coords[j];
        		if (p1 == null)
        			continue;
        		for (int k = i; k < i + this.length - 1; k++) {
        			if (k == j)
        				continue;
        			Point3d p2 = coords[k];
        			if (p2 == null)
        				continue;
        			double dist = p1.distance(p2);
        			features[i] += dist * dist;
        		}
    		}
    		features[i] = features[i]/(this.length * this.length);
    	}
		return new EndToEndDistanceDoubleSequenceFeature(features);
	}
}
