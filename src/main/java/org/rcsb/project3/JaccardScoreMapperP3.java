package org.rcsb.project3;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data list, to
 * a Jaccard Index. It calculates the Jaccard index for multi-sets.
 * 
 * @author  Peter Rose
 */
public class JaccardScoreMapperP3 implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data = null;


	public JaccardScoreMapperP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}

	/**
	 * Returns <PdbId.Chain, Jaccard Index> pairs. This is an extension of the 
	 * Jaccard Index to multi-sets. The multi-sets are represented as vectors,
	 * where each vector element is a feature count.
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,SequenceFeatureInterface<?>> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,SequenceFeatureInterface<?>> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		// get the 2 sequence
		SequenceFeatureInterface<?> v1 = t1._2;
		SequenceFeatureInterface<?> v2 = t2._2;

		float union = 0;
		float intersection = 0;
		int length = v1.length();
		if (v2.length() < length)
			length = v2.length();
		for (int i = 0, n = length; i < n; i++) {
			double v1value = v1.todouble(i);
			double v2value = v2.todouble(i);
			if (v1value < v2value) {
				intersection += v1value;
				union += v2value;
			} else {
				intersection += v2value;
				union += v1value;
			}
		}

		float value = 0;
		if (union > 0) {
		    value = intersection/union;
		}
		
        return new Tuple2<String, Float>(key.toString(), value);
    }
}
