package org.rcsb.rehsDavidM;



	import java.util.HashMap;
	import java.util.List;
	import java.util.Map;

	import javax.vecmath.Point3d;

	import org.apache.spark.broadcast.Broadcast;
	import org.rcsb.project3.*;

	import scala.Tuple2;

	/**
	 * This class maps a pair of chains, specified by two indices into the broadcasted sequences list, to
	 * a Jaccard Index. It calculates the Jaccard index for multi-sets.
	 * 
	 * Order of fragments are getting ignored for Jaccard index calculation.
	 * 
	 * @author  David Mao
	 */
	public class WatermanLevenschteinMapper implements AlignmentAlgorithmInterface {
		private static final long serialVersionUID = 1L;
		private Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences = null;
		double open = 1;
		double extend = 0.1;
	    private int traceback = 0;

		public static void main(String[] args) {
			
		}
		public WatermanLevenschteinMapper() {
		}

		public WatermanLevenschteinMapper(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences) {
			this.sequences = sequences;
		}

		public String getName() {
			return "meetMinIndex";
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
			if (v1 == null || v2 == null){
				return null;
			}
			String key = tuple._1 + "," + tuple._2;
			float value = meetMinIndex(v1, v2);
		
	        return new Tuple2<String, Float>(key, value);
	    }
		
		public Tuple2<String, Float> callLevenshtein(Tuple2<String, String> tuple) throws Exception {
			SequenceFeatureInterface<?> t1 = this.sequences.getValue().get(tuple._1);
			SequenceFeatureInterface<?> t2 = this.sequences.getValue().get(tuple._2);
			
			if (t1 == null || t2 == null) {
				return null;
			}
			
			String key = tuple._1+  "," + tuple._2;
			Float value = (float) LevenshteinDistanceP3.normalizedDistance(t1,  t2);
			
	        return new Tuple2<String, Float>(key.toString(), value);
	    }
		
		public Tuple2<String, Float> callGotoh(Tuple2<String, String> tuple) {
			SequenceFeatureInterface<?> t1 = this.sequences.getValue().get(tuple._1);
			SequenceFeatureInterface<?> t2 = this.sequences.getValue().get(tuple._2);
			
			if (t1 == null || t2 == null) {
				return null;
			}
			
			String key = tuple._1 +"," + tuple._2;
			
			Alignment<?> SWAlignment =getAlignment(t1, t2, open, extend);
			float value = (float) SWAlignment.calculateScore();
			int v1L = t1.length();
			int v2L = t2.length();
			if (v1L > v2L)
				value = value/v2L;
			else
				value = value/v1L;
			// Traceback
			if (traceback > 0)  {
				printTraceback(t1,t2,SWAlignment.getSequence1(),SWAlignment.getSequence2());
			}
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
		public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> coords) {		
		}
		
		private <T> float meetMinIndex(SequenceFeatureInterface<T> s1, SequenceFeatureInterface<T> s2) {
			Map<T, Integer>features1 = calcFeatureCounts(s1);
			Map<T, Integer>features2 = calcFeatureCounts(s2);
			
			return (float) MeetMinIndex.meetMinIndex(features1, features2);
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
		private <T,K> Alignment<T> getAlignment(SequenceFeatureInterface<T> v1,SequenceFeatureInterface<K> v2,double o, double e) {
			return SmithWatermanGotoh.align(v1, (SequenceFeatureInterface<T>)v2, o, e);
		}
		private void printTraceback(SequenceFeatureInterface<?> v1,SequenceFeatureInterface<?> v2,Integer[] v1Order,Integer[] v2Order) {
			Integer c1 = v1Order[0];
			Integer c2 = v2Order[0];
			String commonV1 = "start from " + c1 + "\t";
			String commonV2 = "start from " + c2 + "\t";
			for (int i = 0; i < v1Order.length; i++) {
				c1 = v1Order[i];
				c2 = v2Order[i];
				if (c1 == null) {
					commonV1 += "--- \t";
					commonV2 += "xxx \t";
				} else if (c2 == null) {
					commonV1 += "xxx \t";
					commonV2 += "--- \t";
				} else {
					commonV1 += v1.toString(c1) + " \t";
					commonV2 += v2.toString(c2) + " \t";
				}
			}
			commonV1 += " end at " + c1;
			commonV2 += " end at " + c2;
			System.out.println(commonV1);
			System.out.println(commonV2);
			System.out.println();
		}



}
