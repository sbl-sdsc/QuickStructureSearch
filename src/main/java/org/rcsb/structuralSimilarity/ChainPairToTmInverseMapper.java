package org.rcsb.structuralSimilarity;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class maps a pair of chains with the second chain inverted, specified by two indices into the broadcasted data
 * to a vector of alignment scores
 * 
 * @author  Peter Rose
 */
public class ChainPairToTmInverseMapper implements PairFunction<Tuple2<Integer,Integer>,String,Float[]> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, Point3d[]>>> data = null;


	public ChainPairToTmInverseMapper(Broadcast<List<Tuple2<String,Point3d[]>>> data) {
		this.data = data;
	}

	/**
	 * Returns a chainId pair with the TM scores
	 */
	public Tuple2<String, Float[]> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Point3d[]> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Point3d[]> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		Point3d[] points1 = t1._2;
		Point3d[] points2 = invert(t2._2);

		Float[] scores = TmScorer.getFatCatTmScore(points1, points2);
		
        return new Tuple2<String, Float[]>(key.toString(), scores);
    }
	
	private static Point3d[] invert(Point3d[] points) {
		Point3d[] inverse = new Point3d[points.length];
		for (int i = 0; i < points.length; i++) {
			inverse[i] = points[points.length - i - 1];
		}
		return inverse;
	}
}
