package org.rcsb.projectm;

import java.io.Serializable;

import javax.vecmath.Point3d;

import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.project3.SequenceFingerprint;

/**
 * This class generates a fingerprint (signature) for protein chains based upon
 * the end-to-end distance of fragments
 *
 * @author Peter Rose
 */
public class UsrMomentsFingerprint implements SequenceFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
//	private double MAX_AVE_CA_CA_DISTANCE = 3.3;
//	private double[] moments;
 
    /**
     * Default constructor uses default parameters
     */
    public UsrMomentsFingerprint() {
    }
   
    /**
     * Constructor with all parameters
     * @param length fragment length
     */
 /**   public UsrFingerprint (int length, double binSize) {
        this.length = length;
        this.binSize = binSize;
        this.featureCount = (int)(this.length * this.MAX_AVE_CA_CA_DISTANCE/this.binSize);
	}
*/
    public String getName() {
    	return this.getClass().getSimpleName();
    }

    /**
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
    public UsrFeature getFingerprint(Point3d[] coords) {
    	double[] features = new double[12];
/**
    	double scale = 1/this.binSize;
    	if (coords.length-this.length < 0) {
    		return features;
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
    		int bin = (int)Math.round(scale*first.distance(last));
    		if (bin > this.featureCount-1) {
    			continue;
    		}
    		features[bin]++;
    	}
    	*/
    	features = GenerateMoments.getMoments(coords);
//    	System.out.println(Arrays.toString(features));
		return new UsrFeature(features);
    }
}
