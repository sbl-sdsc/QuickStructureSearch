package org.rcsb.fingerprints;

import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * This class generates a fingerprint (signature) for protein chains based upon
 * the end-to-end distance of fragments
 *
 * @author Peter Rose
 */
public class EndToEndDistanceFingerprint implements GenericFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private double MAX_AVE_CA_CA_DISTANCE = 3.3;
    private int length = 9;
    private double binSize = 2.0;
    private int featureCount = (int)(this.length * this.MAX_AVE_CA_CA_DISTANCE/this.binSize);
 
    /**
     * Default constructor uses default parameters
     */
    public EndToEndDistanceFingerprint() {}
    
    /**
     * Constructor with all parameters
     * @param length fragment length
     */
    public EndToEndDistanceFingerprint (int length, double binSize) {
        this.length = length;
        this.binSize = binSize;
        this.featureCount = (int)(this.length * this.MAX_AVE_CA_CA_DISTANCE/this.binSize);
	}

    public String getName() {
    	return this.getClass().getSimpleName() + "_L" + this.length + "B" + this.binSize;
    }

    /**
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
    public double[] getFingerprint(Point3d[] coords) {
    	double[] features = new double[this.featureCount];

    	double scale = 1/this.binSize;
    	if (coords.length-this.length-1 <= 0) {
    		return features;
    	}
    	for (int i = 0; i < coords.length-this.length-1; i++) {
    		Point3d first = coords[i];
    		Point3d last = coords[i+this.length-1];
    		
    		// skip gaps
    		if (first == null || last  == null) {
    			continue;
    		}

    		// calculate end to end distance of fragment
    		// and bin values
    		int bin = (int)Math.round(scale*first.distance(last));
    		if (bin > this.featureCount-1) {
    			continue;
    		}
    		features[bin]++;
    	}
//    	System.out.println(Arrays.toString(features));
		return features;
    }
}
