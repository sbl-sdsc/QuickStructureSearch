package org.rcsb.projectec;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;
import org.rcsb.project3.AlignmentAlgorithmInterface;
import org.rcsb.project3.SequenceFeatureInterface;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted sequences list, to
 * a Jaccard Index. It calculates the Jaccard index for multi-sets.
 * 
 * Order of fragments are getting ignored for Jaccard index calculation.
 * 
 * @author  Peter Rose, Chris Li
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
		return "NormalizedCompressionDistance";
	}
	
	/**
	 * Returns <PdbId.Chain, Jaccard Index> pairs. This is an extension of the 
	 * Jaccard Index to multi-sets. The multi-sets are represented as vectors,
	 * where each vector element is a feature count.
	 */
	@SuppressWarnings("unchecked")
	public Tuple2<String, Float> call(Tuple2<String, String> tuple) throws Exception {
		SequenceFeatureInterface<Integer> v1 = (SequenceFeatureInterface<Integer>) this.sequences.getValue().get(tuple._1);
		SequenceFeatureInterface<Integer> v2 = (SequenceFeatureInterface<Integer>) this.sequences.getValue().get(tuple._2);

		if (v1 == null || v2 == null) {
			return null;
		}
		
		String key = tuple._1 + "," + tuple._2;
		float value = meetMinIndex(v1, v2);
	
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
	
	private <T> float meetMinIndex(SequenceFeatureInterface<T> s1, SequenceFeatureInterface<T> s2) {
		int[] array1 = toIntArray(s1);
		int[] array2 = toIntArray(s2);
		return (float) NormalizedCompressionDistanceEC.distance(array1, array2);
	}
	
	/**
	 * Counts how many times each feature occurs in the sequence.
	 * @param sequence map where the key represents the feature and the value represents the feature count
	 * @return
	 */
	private <T> int[] toIntArray(SequenceFeatureInterface<T> sequence) {
		int[] array = new int[sequence.length()];
		T[] array2 = sequence.getSequence();
		for (int i = 0; i < array.length; i++) {
			array[i] = (int)array2[i];
		}
		return array;
	}
}
