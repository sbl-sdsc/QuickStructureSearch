package org.rcsb.project10;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data
 * to a vector of alignment scores
 * 
 * @author  Peter Rose
 */
public class ChainPairToTmMapper implements PairFunction<Tuple2<Integer,Integer>,String,Float[]> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, WritableSegment>>> data = null;


	public ChainPairToTmMapper(Broadcast<List<Tuple2<String,WritableSegment>>> data) {
		this.data = data;
	}

	/**
	 * Returns a chainId pair with the TM scores
	 */
	public Tuple2<String, Float[]> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,WritableSegment> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,WritableSegment> t2 = this.data.getValue().get(tuple._2);
		
		String key = t1._1 + "," + t2._1;
//		StringBuilder key = new StringBuilder();
//		key.append(t1._1);
//		key.append(",");
//		key.append(t2._1);
		
		Point3d[] points1 = t1._2.getCoordinates();
		Point3d[] points2 = t2._2.getCoordinates();

		Float[] scores = TmScorer.getFatCatTmScore(points1, points2);
//		Float[] scores = new Float[1];
//		scores[0] = (float)Math.abs(points1.length - points2.length);
		
        return new Tuple2<String, Float[]>(key.toString(), scores);
    }
}