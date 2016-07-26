package org.rcsb.rehsDavidM;
import org.rcsb.project3.*;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * 
 * @author  David Mao
 */
public class SmithWatermanMapperDM implements AlignmentAlgorithmInterface {
	private static final long serialVersionUID = 1L;
	private Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences = null;
	private double start = 1;
	private double extend = 0.1;
	public static void main(String[] args) {
		
	}
	public SmithWatermanMapperDM() {
	}

	public SmithWatermanMapperDM(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}

	public String getName() {
		return "ModifiedSmithWaterman";
	}
	
	/**
	 * 
	 */
	@SuppressWarnings("unchecked")
	public Tuple2<String, Float> call(Tuple2<String, String> tuple) throws Exception {
		SequenceFeatureInterface<Integer> v1 = (SequenceFeatureInterface<Integer>) this.sequences.getValue().get(tuple._1);
		SequenceFeatureInterface<Integer> v2 = (SequenceFeatureInterface<Integer>) this.sequences.getValue().get(tuple._2);
		
		if (v1 == null || v2 == null) {
			return null;
		}

		String key = tuple._1 + "," + tuple._2;
		float value = hashingJaccardIndex(v1, v2);
        return new Tuple2<String, Float>(key, value);
        
    }

	@Override
	public void setSequence(Broadcast<Map<String, SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}
//
	/**
	 * Not used in this algorithm
	 */
	@Override
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> coords) {		
	}
	
	/**
	 * Hashes two SequenceFeatureInterfaces to int arrays. 
	 * Uses generalized Jaccard distance to compute distance.
	 * @param s1
	 * @param s2
	 * @return
	 */
	private <T> float hashingJaccardIndex(SequenceFeatureInterface<T> s1, SequenceFeatureInterface<T> s2) {
//		Map<T, Integer>features1 = calcFeatureCounts(s1);
//		Map<T, Integer>features2 = calcFeatureCounts(s2);
//		
//		return (float) HashingJaccardIndex.HashingJaccardIndex(features1, features2);
		int [] seq1 = new int[s1.length()];
		int [] seq2 = new int [s2.length()];
		for (int i = 0; i < seq1.length; i++){
			seq1[i] = (int) s1.get(i);
		}
		for (int i = 0; i < seq2.length; i++){
			seq2[i] = (int) s2.get(i);
		}
		
		return (float)SmithWatermanDM.distance(seq1, seq2,start,extend);
		
		
	}
	
	/**
	 * Counts how many times each feature occurs in the sequence.
	 * @param sequence map where the key represents the feature and the value represents the feature count
	 * @return
	 */
//	private <T> Map<T, Integer> calcFeatureCounts(SequenceFeatureInterface<T> sequence) {
//		Map<T, Integer> featureCounts = new HashMap<T, Integer>(sequence.length());
//		for (T key: sequence.getSequence()) {
//			Integer count = featureCounts.getOrDefault(key, 0);
//			featureCounts.put(key, count+1);
//		}
//		return featureCounts;
//	}
	
	
}
