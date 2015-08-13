package org.rcsb.project2;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import scala.Tuple2;

public class SecondaryStructFeature implements StructureCompare {
	protected Point3d[] pts;
	protected Tuple2<int[], int[]> features;
	protected SecondaryStructProjection[] projections;
	protected Vector3d normP, normX;
	protected Point3d normC;

	@Override
	public SecondaryStructProjection getNormProjection(byte b) {
		Vector3d x = new Vector3d(normX);
		Vector3d p = new Vector3d(normP);
		// System.out.printf("s1: %d,\ts2: %d" + System.lineSeparator(), ((b % 2) << 1) - 1, ((b >> 1) << 1) - 1);
		x.scale(((b % 2) << 1) - 1);
		p.scale(((b >> 1) << 1) - 1);
		return SecondaryStructTools.projectOnto(pts, features, p, new Point3d(normC), x);
	}

	@Override
	public SecondaryStructProjection getProjection(int i) {
		return projections[i] == null ? projections[i] = SecondaryStructTools.project(pts, features, i)
				: projections[i];
	}

	@Override
	public int length() {
		return features._1.length;
	}

	public int lengthOf(int i) {
		return features._2[i] - features._1[i];
	}

	public Tuple2<int[], int[]> getFeatures() {
		return features;
	}

	public int getNumPoints() {
		int o = 0;
		for (int i = 0; i < length(); i++)
			o += features._2[i] - features._1[i];
		return o;
	}
}
