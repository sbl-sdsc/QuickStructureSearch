package org.rcsb.structuralSimilarity;

import java.util.Arrays;
import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data list, to
 * a Containment Score. It calculates the ContainmentScore for multi-sets.
 * 
 * @author  Peter Rose
 */
public class FeatureVectorToCosineScoreMapper implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,Vector>>> data = null;


	public FeatureVectorToCosineScoreMapper(Broadcast<List<Tuple2<String,Vector>>> data) {
		this.data = data;
	}

	/**
	 * Returns <PdbId.Chain, Containment Score> pairs. This is an extension of the 
	 * Containment socre to multi-sets. The multi-sets are represented as vectors,
	 * where each vector element is a feature count.
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Vector> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Vector> t2 = this.data.getValue().get(tuple._2);

		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);

		double[] v1 = t1._2.toArray();
		double[] v2 = t2._2.toArray();
		
//		System.out.println("1: " + Arrays.toString(v1));
//		System.out.println("2: " + Arrays.toString(v2));
		
		v1 = normalize(v1);
		v2 = normalize(v2);

		double cos = 0;
	    for (int i = 0; i < v1.length; i++) {
	    	cos += v1[i] * v2[i];
	    }
//	    System.out.println("cos: " + cos);
	    
		return new Tuple2<String, Float>(key.toString(), (float)cos);
	}
	
	private static double[] normalize(double[] vector) {
		double[] vector2 = new double[vector.length];
		double sum = 0;

		for (double v: vector2) {
			sum += v*v;
		}
		for (int i = 0; i < vector.length; i++) {
			vector2[i] = vector2[i]/sum;
		}

		return vector2;
	}
}
