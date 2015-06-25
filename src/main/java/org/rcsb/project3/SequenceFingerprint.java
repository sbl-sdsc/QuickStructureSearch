package org.rcsb.project3;

import javax.vecmath.Point3d;

/**
 * This interface is used for fingerprint that return as a SequenceFeature
 * 
 * @author Chris Li
 */
public interface SequenceFingerprint {
	public SequenceFeatureInterface<?> getFingerprint(Point3d[] coordinates);
	public String getName();
}
