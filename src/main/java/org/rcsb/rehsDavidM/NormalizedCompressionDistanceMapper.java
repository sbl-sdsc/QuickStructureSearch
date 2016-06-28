package org.rcsb.rehsDavidM;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;
import org.rcsb.project3.*;

import scala.Tuple2;

/**
 * 
 * @author  David Mao
 */
public class NormalizedCompressionDistanceMapper implements AlignmentAlgorithmInterface {
	private static final long serialVersionUID = 1L;
	private Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences = null;

	public static void main(String[] args) {
		
	}
	public NormalizedCompressionDistanceMapper() {
	}

	public NormalizedCompressionDistanceMapper(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}

	public String getName() {
		return "meetMinIndex";
	}
	
	/**
	 * Returns <
	 * Returns <PdbId.Chain, Jaccard Index> pairs. This is an extension of the 
	 * Jaccard Index to multi-sets. The multi-sets are represented as vectors,
	 * where each vector element is a feature count.
	 */
	@SuppressWarnings("unchecked")
	public Tuple2<String, Float> call(Tuple2<String, String> tuple) throws Exception {
		SequenceFeatureInterface<Integer> v1 = (SequenceFeatureInterface<Integer>) this.sequences.getValue().get(tuple._1);
		SequenceFeatureInterface<Integer> v2 = (SequenceFeatureInterface<Integer>) this.sequences.getValue().get(tuple._2);
		String key = tuple._1 + "," + tuple._2;
		float value = (float) distance(v1, v2);
	
        return new Tuple2<String, Float>(key, value);
    }

	@Override
	public void setSequence(Broadcast<Map<String, SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}

	/**
	 * Not used in this algorithm
	 */
	@Override
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> coords) {		
	}
	
	private <T> float distance(SequenceFeatureInterface<T> s1, SequenceFeatureInterface<T> s2) {
		int[] w1 = new int[s1.length()];
		for (int i = 0; i < s1.length(); i++){
			w1[i] = (int) s1.get(i);
		}
		int[] w2 = new int[s2.length()];
		for (int i = 0; i < s2.length(); i++){
			w2[i] = (int) s2.get(i);
		}
		return (float) NormalizedCompressionDistance.distance(w1,w2);
	}
	
	/**
	 * Counts how many times each feature occurs in the sequence.
	 * @param sequence map where the key represents the feature and the value represents the feature count
	 * @return
	 */
	private <T> Map<T, Integer> calcFeatureCounts(SequenceFeatureInterface<T> sequence) {
		Map<T, Integer> featureCounts = new HashMap<T, Integer>(sequence.length());
		for (T key: sequence.getSequence()) {
			Integer count = featureCounts.getOrDefault(key, 0);
			featureCounts.put(key, count+1);
		}
		return featureCounts;
	}
}
