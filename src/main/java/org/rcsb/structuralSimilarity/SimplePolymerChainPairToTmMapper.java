package org.rcsb.structuralSimilarity;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.hadoop.io.SimplePolymerChain;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data
 * to a vector of alignment scores
 * 
 * @author  Peter Rose
 */
public class SimplePolymerChainPairToTmMapper implements PairFunction<Tuple2<Integer,Integer>,String,Float[]> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, SimplePolymerChain>>> data = null;


	public SimplePolymerChainPairToTmMapper(Broadcast<List<Tuple2<String, SimplePolymerChain>>> chainsBc) {
		this.data = chainsBc;
	}

	/**
	 * Returns a chainId pair with the TM scores
	 */
	public Tuple2<String, Float[]> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,SimplePolymerChain> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,SimplePolymerChain> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		Point3d[] points1 = t1._2.getCoordinates();
		Point3d[] points2 = t2._2.getCoordinates();

		Float[] scores = TmScorer.getFatCatTmScore(points1, points2);
		
        return new Tuple2<String, Float[]>(key.toString(), scores);
    }
}