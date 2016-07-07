package org.rcsb.projectva;
import org.rcsb.project3.*;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class compares and scores a pair of chains' fingerprint sequences with the Levenshtein Distance algorithm.
 * 
 * @author  Peter Rose, Chris Li
 */
public class LevenshteinMapperva implements AlignmentAlgorithmInterface {
	private static final long serialVersionUID = 1L;
	private Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences = null;

	public LevenshteinMapperva() {
	}

	public LevenshteinMapperva(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}
	
	public String getName() {
		return "LevenshteinIndex";
	}
	
	/**
	 * Returns <PdbId.Chain, Levenshtein distance> pairs. 
	 */
	public Tuple2<String, Float> call(Tuple2<String, String> tuple) throws Exception {
		SequenceFeatureInterface<?> t1 = this.sequences.getValue().get(tuple._1);
		SequenceFeatureInterface<?> t2 = this.sequences.getValue().get(tuple._2);
		
		if (t1 == null || t2 == null) {
			return null;
		}
		
		String key = tuple._1+  "," + tuple._2;
		Float value = (float) LevenshteinDistanceP3.normalizedDistance(t1,  t2);
		
        return new Tuple2<String, Float>(key.toString(), value);
    }

	@Override
	public void setSequence(Broadcast<Map<String, SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}

	/**
	 * Not used in this algorithm
	 */
	@Override
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> sequence) {
	}
}
