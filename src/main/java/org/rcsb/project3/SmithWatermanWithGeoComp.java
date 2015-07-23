package org.rcsb.project3;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the longest local common subsequence over the length of the chains
 * using SmithWaterman algorithm and Gotoh's improvement
 * 
 * @author Chris Li
 */
public class SmithWatermanWithGeoComp implements AlignmentAlgorithmInterface {

	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data = null;
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

    public SmithWatermanWithGeoComp() {
	}
    
    public SmithWatermanWithGeoComp(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
    
    /***
     * Constructor with setting options
     * @param data
     * @param open
     * @param extend
     */
    public SmithWatermanWithGeoComp(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data, double open, double extend) {
		this.data = data;
		this.open = open;
		this.extend = extend;
	}
    
    /**
     * Constructor with traceback option
     * @param data
     * @param traceback
     */
	public SmithWatermanWithGeoComp(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data, int traceback) {
		this.data = data;
		this.traceback = traceback;
	}
	
	@Override
	public void setSequence(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
	
	@Override
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) {
		Tuple2<String,SequenceFeatureInterface<?>> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,SequenceFeatureInterface<?>> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		SequenceFeatureInterface<?> v1 = t1._2;
		SequenceFeatureInterface<?> v2 = t2._2;
		
		Alignment<?> SWAlignment = getAlignment(v1, v2, open, extend);
		
		Point3d[] c1 = v1.getCoords();
		Point3d[] c2 = v2.getCoords();

		Integer[] v1Order = SWAlignment.getSequence1();
		Integer[] v2Order = SWAlignment.getSequence2();
		if (v1Order.length > 0) {
			int min = Integer.MAX_VALUE;
			int minIndex = 0;
			for (int i = 0; i < v1Order.length; i++) {
				Integer a1 = v1Order[i];
				Integer a2 = v2Order[i];
				if (a1 != null && a2 != null) {
					int unAlign = 0;
					for (int j =0; j < v1Order.length; j++) {
						Integer b1 = v1Order[j];
						Integer b2 = v2Order[j];
						if (b1 != null && b2 != null) {
							double dist1 = c1[a1].distance(c1[b1]);
							double dist2 = c2[a2].distance(c2[b2]);
							if (dist1 - dist2 < -8 || dist1 - dist2 > 8) {
								unAlign++;
							}
						}
					}
					if (unAlign < min) {
						min = unAlign;
						minIndex = i;
					}
				}
			}
			Integer min1 = v1Order[minIndex];
			Integer min2 = v2Order[minIndex];
			for (int i =0; i < v1Order.length; i++) {
				Integer a1 = v1Order[i];
				Integer a2 = v2Order[i];
				if (a1 != null && a2 != null) {
					double dist1 = c1[min1].distance(c1[a1]);
					double dist2 = c2[min2].distance(c2[a2]);
					if (dist1 - dist2 < -8 || dist1 - dist2 > 8) {
						if ((i != 0 && v1Order[i-1] == null) || (i != v1Order.length-1 && v1Order[i+1] == null))
							v1Order[i] = null;
						else
							v2Order[i] = null;
	 				}
				}
			}
		}
		
		SWAlignment.setSequence1(v1Order);
		SWAlignment.setSequence2(v2Order);
		
		float value = (float) SWAlignment.calculateScore();
		int v1L = v1.length();
		int v2L = v2.length();
		if (v1L > v2L)
			value = value/v2L;
		else
			value = value/v1L;
		// Traceback
		if (traceback > 0)  {
			printTraceback(v1,v2,SWAlignment.getSequence1(),SWAlignment.getSequence2());
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
}
