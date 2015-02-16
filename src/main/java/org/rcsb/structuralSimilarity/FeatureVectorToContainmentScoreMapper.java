package org.rcsb.structuralSimilarity;

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
public class FeatureVectorToContainmentScoreMapper implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,Vector>>> data = null;


	public FeatureVectorToContainmentScoreMapper(Broadcast<List<Tuple2<String,Vector>>> data) {
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

		float set1 = 0;
		float set2 = 0;
		float intersection = 0;
		for (int i = 0, n = v1.length; i < n; i++) {
			set1 += v1[i];
			set2 += v2[i];
			if (v1[i] < v2[i]) {
				intersection += v1[i];
			} else {
				intersection += v2[i];
			}
		}

		float value = 0;
		if (set1 > 0 && set2 > 0) {
			if (set1 > set2) {
				value = intersection/set2;
			} else {
				value = intersection/set1;
			}
		}

		return new Tuple2<String, Float>(key.toString(), value);
	}
}
