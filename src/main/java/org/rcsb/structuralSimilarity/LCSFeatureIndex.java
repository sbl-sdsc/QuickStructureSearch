package org.rcsb.structuralSimilarity;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the length longest common subsequence over the length of the chains
 */
public class LCSFeatureIndex implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	// the amount to consider two value different
	private double diff = 3.14/72;
	private Broadcast<List<Tuple2<String,Vector>>> data = null;

	/**
     * Constructor with default setting
     */
	public LCSFeatureIndex(Broadcast<List<Tuple2<String,Vector>>> data) {
		this.data = data;
	}
	
	/**
     * Constructor with parameter
     */
	public LCSFeatureIndex(Broadcast<List<Tuple2<String,Vector>>> data, double diff) {
		this.data = data;
		this.diff = diff;
	}

	/**
	 * Returns <PdbId.Chain, LCS score> pairs.
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Vector> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Vector> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		// get the 2 sequence
		double[] v1 = t1._2.toArray();
		double[] v2 = t2._2.toArray();
		
		// The longest common substring
		float LCSres = (float)LCS(v1,v2,diff);
		float value = LCSres/(float)getLength(v1,v2);
        return new Tuple2<String, Float>(key.toString(), value);
    }
	
	/**
	 * Get the length of the longest common substring
	 */
	private static int LCS(double[] v1, double[] v2, double diff) {
		// the length of LCS at each pair of values
		int c[][] = new int[v1.length][v2.length];
		// the first row
		for (int i = 0; i < v1.length; i++) {
			if (similar(v1[i],v2[0],diff))
				c[i][0] = 1;
			else 
				c[i][0] = 0;
		}
		// the first column
		for (int j = 0; j < v2.length; j++) {
			if (similar(v1[0],v2[j],diff)) {
				c[0][j] = 1;
			}
			else 
				c[0][j] = 0;
		}
		// the rest
		for (int i = 1; i < v1.length; i++) {
			for (int j = 1; j < v2.length; j++) {
				if (similar(v1[i],v2[j],diff)) {
					c[i][j] = c[i-1][j-1] + 1;
				}
				else if (c[i-1][j] >= c[i][j-1]) {
					c[i][j] = c[i-1][j];
				}
				else {
					c[i][j] = c[i][j-1];
				}
			}
		}
		// the last pair represent the global LCS
		return c[v1.length-1][v2.length-1];
	}
	
	/**
	 * Compare the value
	 */
	private static boolean similar(double x, double y, double diff) {
		// if the value is NaN, it is a gap
		if (Double.isNaN(x) || Double.isNaN(y)){
			return false;
		}
		// if the difference is smaller than the setting diff, the pairs are similar 
		if (((x+diff) > y) && ((x-diff) < y))
			return true;
		else return false;
	}
	
	/**
	 * Return the length of the shorter sequence
	 */
	private static int getLength(double[] v1, double[] v2) {
		int v1l = v1.length;
		int v2l = v2.length;
		// gap is not counted
		for (int i = 0; i < v1.length; i++) {
			if (Double.isNaN(v1[i]))
				v1l--;
		}
		for (int i = 0; i < v2.length; i++) {
			if (Double.isNaN(v2[i]))
				v2l--;
		}
		if (v1l < v2l)
			return v1l;
		else return v2l;
	}
}
