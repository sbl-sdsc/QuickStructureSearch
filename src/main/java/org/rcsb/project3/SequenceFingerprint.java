package org.rcsb.project3;

import javax.vecmath.Point3d;

/**
 * This interface is used for fingerprint that return as a SequenceFeature
 * 
 * @author Chris Li
 */
public interface SequenceFingerprint {
	/**
	 * transfer the Point3d coordinates into a SequenceFeature fingerprint
	 * @param coordinates
	 * @return
	 */
	public SequenceFeatureInterface<?> getFingerprint(Point3d[] coordinates);
	/**
	 * get the name of the fingerprint
	 * @return
	 */
	public String getName();
}
