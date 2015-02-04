package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

public class ChainQualityMapper implements PairFunction<Tuple2<String,Point3d[]>,String, Integer> {
	private static final long serialVersionUID = 1L;

	@Override
	public Tuple2<String, Integer> call(Tuple2<String, Point3d[]> t)
			throws Exception {
        Point3d[] points = t._2;
		
		int length = 0;
		
		// TODO exclude his tags
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				length++;
			}
		}
		
		return new Tuple2<String, Integer>(t._1, length);
	}

}
