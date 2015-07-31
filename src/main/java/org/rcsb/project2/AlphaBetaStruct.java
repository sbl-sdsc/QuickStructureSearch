package org.rcsb.project2;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import scala.Tuple2;

public class AlphaBetaStruct extends SecondaryStructFeature {
	public AlphaBetaStruct(Tuple2<int[], int[]> features, Vector3d normP, Vector3d normX, Point3d normC, Point3d[] pts) {
		this.features = features;
		this.normP = normP;
		this.normX = normX;
		this.normC = normC;
		this.pts = pts;
		this.projections = new SecondaryStructProjection[length()];
	}
}
