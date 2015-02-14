package org.rcsb.structuralSimilarity;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class filters chain pairs by difference in chain length.
 * 
 * @author  Peter Rose
 */
public class ChainPairLengthFilter implements Function<Tuple2<Integer,Integer>, Boolean> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, Point3d[]>>> data = null;
	private double coverage;

/**
 * @param data Broadcast list of <ChainId, Calpha coordinate array> pairs
 * @param coverage fraction of the short chain that covers the long chain
 */
	public ChainPairLengthFilter(Broadcast<List<Tuple2<String,Point3d[]>>> data, double coverage) {
		this.data = data;
		this.coverage = coverage;
	}

	/**
	 * Returns true if the length of the shorter chain is at least a given fraction of the length of the longer chain
	 */
	public Boolean call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Point3d[]> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Point3d[]> t2 = this.data.getValue().get(tuple._2);
		
		int minLen = Math.min(t1._2.length, t2._2.length);
		int maxLen = Math.max(t1._2.length, t2._2.length);
		
		return (minLen/(double)maxLen >= coverage);
    }
}
