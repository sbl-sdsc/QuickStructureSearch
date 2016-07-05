package org.rcsb.projectm.Demo;

import java.util.List;

import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;
import org.rcsb.project3.*;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the longest local common subsequence over the length of the chains
 * using SmithWaterman algorithm and Gotoh's improvement
 * 
 * @author Chris Li
 */
public class SmithWatermanGotohMapperP3 implements AlignmentAlgorithmInterface {

	private static final long serialVersionUID = 1L;
	private Broadcast<Map<String,SequenceFeatureInterface<?>> >sequences = null;
    // print traceback if it is greater than 0
    private int traceback = 0;
    /* With different open and extend penalty, this class could function the same as LCS or SmithWaterman
     * LCS: open = extend = 0;
     * SmithWaterman = open = extend = 1;
     */
    // open gap penalty
    private double open = 1;
    // extend gap penalty
    private double extend = 0.1;

    public SmithWatermanGotohMapperP3() {
	}
    
    public SmithWatermanGotohMapperP3(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}
    
	public String getName() {
		return "SmithWatermanGotoh";
	}
    
    /***
     * Constructor with setting options
     * @param sequences
     * @param open
     * @param extend
     */
    public SmithWatermanGotohMapperP3(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences, double open, double extend) {
		this.sequences = sequences;
		this.open = open;
		this.extend = extend;
	}
    
    /**
     * Constructor with traceback option
     * @param sequences
     * @param traceback
     */
	public SmithWatermanGotohMapperP3(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences, int traceback) {
		this.sequences = sequences;
		this.traceback = traceback;
	}
	
	@Override
	public void setSequence(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences) {
		this.sequences = sequences;
	}
	
	@Override
	public Tuple2<String, Float> call(Tuple2<String, String> tuple) {
		SequenceFeatureInterface<?> t1 = this.sequences.getValue().get(tuple._1);
		SequenceFeatureInterface<?> t2 = this.sequences.getValue().get(tuple._2);
		
		if (t1 == null || t2 == null) {
			return null;
		}
		
		String key = tuple._1 +"," + tuple._2;
		
		Alignment<?> SWAlignment = getAlignment(t1, t2, open, extend);
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
	
	/**
	 * Get alignment for the two sequence. Object class casting.
	 * @param v1
	 * @param v2
	 * @param o
	 * @param e
	 * @return
	 */
	@SuppressWarnings("unchecked")
	private <T,K> Alignment<T> getAlignment(SequenceFeatureInterface<T> v1,SequenceFeatureInterface<K> v2,double o, double e) {
		return SmithWatermanGotoh.align(v1, (SequenceFeatureInterface<T>)v2, o, e);
	}
	
	/**
	 * Print the SmithWatermanGotoh traceback
	 * @param v1
	 * @param v2
	 * @param b
	 */
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

	/**
	 * Not used in this algorithm
	 */
	@Override
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> sequence) {
	}
}
