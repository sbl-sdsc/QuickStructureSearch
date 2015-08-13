package org.rcsb.project3;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class compares and scores a pair of chains' fingerprint sequences with the Levenshtein Distance algorithm.
 * 
 * @author  Peter Rose, Chris Li
 */
public class LevenshteinMapperP3 implements AlignmentAlgorithmInterface {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> sequences = null;

	public LevenshteinMapperP3() {
	}

	public LevenshteinMapperP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> sequences) {
		this.sequences = sequences;
	}

	/**
	 * Returns <PdbId.Chain, Levenshtein distance> pairs. 
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,SequenceFeatureInterface<?>> t1 = this.sequences.getValue().get(tuple._1);
		Tuple2<String,SequenceFeatureInterface<?>> t2 = this.sequences.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		SequenceFeatureInterface<?> v1 = t1._2;
		SequenceFeatureInterface<?> v2 = t2._2;

		Float value = (float) LevenshteinDistanceP3.normalizedDistance(v1,  v2);
		
        return new Tuple2<String, Float>(key.toString(), value);
    }

	@Override
	public void setSequence(Broadcast<List<Tuple2<String, SequenceFeatureInterface<?>>>> sequences) {
		this.sequences = sequences;
	}

	/**
	 * Not used in this algorithm
	 */
	@Override
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> sequence) {
	}
}
