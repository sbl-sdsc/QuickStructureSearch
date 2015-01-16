package org.rcsb.fingerprints;

import javax.vecmath.Point3d;

public interface GenericFingerprint {
	public double[] getFingerprint(Point3d[] coordinates);
	public String getName();
}
