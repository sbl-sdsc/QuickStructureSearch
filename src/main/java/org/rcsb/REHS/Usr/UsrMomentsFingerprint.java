package org.rcsb.REHS.Usr;

import java.io.Serializable;

import javax.vecmath.Point3d;

import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.project3.SequenceFingerprint;

/**
 * This class generates a fingerprint (signature) for protein chains based upon
 * its USR moments.
 *
 * @author Peter Rose, Michael Wang
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
    	features = GenerateMoments.getMoments(coords);
		return new UsrFeature(features);
    }
}
