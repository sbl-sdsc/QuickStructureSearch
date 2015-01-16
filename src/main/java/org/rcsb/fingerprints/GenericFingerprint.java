package org.rcsb.fingerprints;

import javax.vecmath.Point3d;

public interface GenericFingerprint {
	/**
	 * Returns a feature vector that represents feature counts, i.e.,
	 * each element of the vector is the number of time a feature
	 * occurs in the macromolecule chain, represent by it 3D coordinates.
	 * @param coordinates
	 * @return feature vector
	 */
	public double[] getFingerprint(Point3d[] coordinates);
	public String getName();
}
