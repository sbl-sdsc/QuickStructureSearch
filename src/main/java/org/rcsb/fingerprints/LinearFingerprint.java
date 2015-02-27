package org.rcsb.fingerprints;

import javax.vecmath.Point3d;

public interface LinearFingerprint {
	/**
	 * Returns a feature vector that represents the features
	 * of a protein chain in the order of the chain.
	 * @param coordinates
	 * @return feature vector
	 */
	public int[] getLinearFingerprint(Point3d[] coordinates);
}
