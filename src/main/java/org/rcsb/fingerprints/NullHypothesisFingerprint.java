package org.rcsb.fingerprints;

import java.io.Serializable;

import javax.vecmath.Point3d;

/**
 * This class generates a fingerprint (signature) for protein chains for null hypothesis testing.
 * The null hypothesis H0 is: all protein fragments of equal length are considered equal for 
 * distinguishing similar from dissimilar proteins. To reject the null hypothesis, the statistical 
 * significance (i.e. 95% or 99% likelihood) of an alternative fingerprint method must be demonstrated.
 * See for example: 
 * https://explorable.com/null-hypothesis
 * http://www.null-hypothesis.co.uk/science/item/what_is_a_null_hypothesis/
 *
 * @author Peter Rose
 */
public class NullHypothesisFingerprint implements GenericFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private int length = 9;
  
    /**
     * Default constructor uses default parameters
     */
    public NullHypothesisFingerprint() {}
    
    /**
     * Constructor with all parameters
     * @param length fragment length
     */
    public NullHypothesisFingerprint (int length) {
        this.length = length;
	}
    
    public String getName() {
    	return this.getClass().getSimpleName() + "_L" + length;
    }
    
    /**
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
    public double[] getFingerprint(Point3d[] coords) {
        double[] fingerPrint = new double[1];
		if (coords.length-length-1 <= 0) {
			return fingerPrint;
		}
		for (int i = 0; i < coords.length-length-1; i++) {
			// skip gaps
    		if (hasGaps(coords, i)) {
    			continue;
    		}
    		// there is only one unique fragment
            fingerPrint[0]++;
		}
		return fingerPrint;
	}
	/**
	 * Returns true if there is a gap between the C alpha atoms
	 * within a fragment. Note, the first position is not checked, 
	 * since the gap info is for a gap between residue i and i + 1.
	 * @param gaps true if there is a gap in the fragment
	 * @param index start residue of fragment
	 * @return
	 */
	private boolean hasGaps(Point3d[] coords, int index) {
		for (int i = index; i < index+this.length; i++) {
			if (coords[i] == null) {
				return true;
			}
		}
		return false;
	}	
}
