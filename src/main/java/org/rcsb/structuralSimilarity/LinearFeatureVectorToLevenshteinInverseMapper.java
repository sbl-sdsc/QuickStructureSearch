package org.rcsb.structuralSimilarity;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class maps a pair of chains, one of them in inverse order, specified by two indices into the broadcasted data list, to
 * a Jaccard Index. It calculates the Jaccard index for multi-sets.
 * 
 * @author  Peter Rose
 */
public class LinearFeatureVectorToLevenshteinInverseMapper implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,Vector>>> data = null;


	public LinearFeatureVectorToLevenshteinInverseMapper(Broadcast<List<Tuple2<String,Vector>>> data) {
		this.data = data;
	}

	/**
	 * Returns <PdbId.Chain, Levenshtein distance> pairs. 
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Vector> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Vector> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		int[] v1 = toIntVector(t1._2.toArray());
		int[] v2 = toInverseIntVector(t2._2.toArray());

		Float value = (float) LevenshteinDistance.normalizedDistance(v1,  v2);
		
        return new Tuple2<String, Float>(key.toString(), value);
    }
	
	private int[] toIntVector(double[] vector) {
		int[] intVector = new int[vector.length];
		for (int i = 0; i < vector.length; i++) {
			intVector[i] = (int)Math.round(vector[i]);
		}
		return intVector;
	}
	
	private int[] toInverseIntVector(double[] vector) {
		int[] intVector = new int[vector.length];
		for (int i = 0; i < vector.length; i++) {
			intVector[i] = (int)Math.round(vector[vector.length - i - 1]);
		}
		return intVector;
	}
	
}
