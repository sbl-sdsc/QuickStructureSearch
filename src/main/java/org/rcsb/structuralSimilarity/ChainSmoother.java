package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

public interface ChainSmoother {
	Point3d[] getSmoothedPoints(Point3d[] points);
}
