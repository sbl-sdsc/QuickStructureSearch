package org.rcsb.project10;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data
 * to a vector of alignment scores
 * 
 * @author  Peter Rose
 */
public class SegmentPairMapperToTmMapper implements PairFunction<Tuple2<Tuple2<String, WritableSegment>, Tuple2<String, WritableSegment>>, String,Float[]> {
	private static final long serialVersionUID = 1L;
	

	/**
	 * Returns a chainId pair with the TM scores
	 */
	public Tuple2<String, Float[]> call(Tuple2<Tuple2<String, WritableSegment>, Tuple2<String, WritableSegment>> t) throws Exception {
		StringBuilder key = new StringBuilder();
		key.append(t._1._1);
		key.append(",");
		key.append(t._2._1);
		
		Point3d[] points1 = t._1._2.getCoordinates();
		Point3d[] points2 = t._2._2.getCoordinates();

		Float[] scores = TmScorer.getFatCatTmScore(points1, points2);
//		Float[] scores = new Float[1];
//		scores[0] = (float)Math.abs(points1.length - points2.length);
		
        return new Tuple2<String, Float[]>(key.toString(), scores);
    }
}