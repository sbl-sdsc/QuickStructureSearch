package org.rcsb.structuralSimilarity;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class maps a pair of chain to the a score from combination of SmithWaterman and LCS
 */
public class FeatureCombination implements PairFunction<Tuple2<Integer,Integer>,String,Float> {

	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,Vector>>> data = null;
	// the amount to consider two value different
	private double diff = 3.14/72;
	// score of match
	private double match = 1;
	// score of mismatch
    private double mismatch = -1;
    // score of indel
    private double indel = -0.1;
    // score of gap
    private double gap = -1;
    // length of v1
    private int v1length;
    // length of v2
    private int v2length;
    // score matrix
    private double [][] score;

    /**
     * Constructor with default setting
     */
	public FeatureCombination(Broadcast<List<Tuple2<String,Vector>>> data) {
		this.data = data;
	}
	
	/**
     * Constructor with parameters
     */
	public FeatureCombination(Broadcast<List<Tuple2<String,Vector>>> data, double diff, 
			double match, double mismatch, double indel, double gap) {
		this.data = data;
		this.match = match;
		this.mismatch = mismatch;
		this.indel = indel;
		this.gap = gap;
	}
	
	/**
	 * Returns pairs of PdbId.Chain and a combination of SmithWaterman and LCS.
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) {
		Tuple2<String,Vector> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Vector> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		// get the 2 sequence
		double[] v1 = t1._2.toArray();
		double[] v2 = t2._2.toArray();
		
		v1length = v1.length;
		v2length = v2.length;
		
		// the length of the shorter chain
		int length = getLength(v1,v2);
		// the SmithWaterman score
		float SmithWatermanres = SmithWatermanScore(v1, v2);
		// the length of LCS
		float LCSres = (float)LCS(v1,v2,diff);
		// Calculate the result value based on the SmithWaterman and LCS
		float value;
		if (LCSres < SmithWatermanres)
			value = LCSres/length;
		else 
			value = SmithWatermanres/length;
		return new Tuple2<String, Float>(key.toString(), value);
	}
	
	/**
	 * Get the length of the longest common substring
	 */
	private int LCS(double[] v1, double[] v2, double diff) {
		// the length of LCS at each pair of values
		int c[][] = new int[v1length][v2length];
		// the first row
		for (int i = 0; i < v1length; i++) {
			if (similar2(v1[i],v2[0],diff))
				c[i][0] = 1;
			else 
				c[i][0] = 0;
		}
		// the first column
		for (int j = 0; j < v2length; j++) {
			if (similar2(v1[0],v2[j],diff)) {
				c[0][j] = 1;
			}
			else 
				c[0][j] = 0;
		}
		// the rest
		for (int i = 1; i < v1length; i++) {
			for (int j = 1; j < v2length; j++) {
				if (similar2(v1[i],v2[j],diff)) {
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
		return c[v1length-1][v2length-1];
	}
	
	/**
	 * Get the score of the pair using SmithWaterman algorithm
	 */
	private float SmithWatermanScore(double[] v1, double[] v2) {
		score = new double[v1length+1][v2length+1];
		score[0][0] = 0;
		// the first row
		for (int i = 1; i <= v1length; i++) {
			score[i][0] = 0;
		}
		// the first column
		for (int j = 1; j <= v2length; j++) {
			score[0][j] = 0;
		}
		// the rest
		for (int i = 1; i <= v1length; i++) {
		    for (int j = 1; j <= v2length; j++) {
		    	// max of the diagonal, up and left
		    	double diagScore = score[i - 1][j - 1] + similarity1(v1[i-1], v2[j-1]);
		    	double upScore = score[i][j - 1] + indel;
		    	double leftScore = score[i - 1][j] + indel;
		    	score[i][j] = maximum(diagScore, upScore,leftScore, 0);
		    }
		}
		return (float) getMaxScore();
	}
	
	/**
	 * Compare the value, return score of different situation 
	 */
	private double similarity1(double x, double y) {
		// if the value is NaN, return gap score
		if (Double.isNaN(x) || Double.isNaN(y)){
			return gap;
		}
		// if the difference is smaller than the setting diff, return match score
		if (((x+diff) > y) && ((x-diff) < y))
			return match;
		else return mismatch;
    }
	
	/**
	 * Compare the value, return true for similar, false for different 
	 */
	private boolean similar2(double x, double y, double diff) {
		if (Double.isNaN(x) || Double.isNaN(y)){
			return false;
		}
		if (((x+diff) > y) && ((x-diff) < y))
			return true;
		else return false;
	}
	
	/**
     * Get the maximum value in the score matrix.
     */
    private double getMaxScore() {
		double maxScore = 0;
		// skip the first row and column
		for (int i = 1; i <= v1length; i++) {
		    for (int j = 1; j <= v2length; j++) {
				if (score[i][j] > maxScore) {
				    maxScore = score[i][j];
				}
		    }
		}
		return maxScore;
    }
    
    /**
     * Maximum of the three
     */
	private double maximum(double a, double b, double c, double d) {
		if (a > b) {
			if (a > c) {
				return a > d ? a : d;
			} else {
				return c > d ? c : d;
			}
		} else if (b > c) {
			return b > d ? b : d;
		} else {
			return c > d ? c : d;
		}
	}
	
	/**
	 * Return the length of the shorter sequence
	 */
	private int getLength(double[] v1, double[] v2) {
		int v1l = v1length;
		int v2l = v2length;
		// gap is not counted
		for (int i = 0; i < v1length; i++) {
			if (Double.isNaN(v1[i]))
				v1l--;
		}
		for (int i = 0; i < v2length; i++) {
			if (Double.isNaN(v2[i]))
				v2l--;
		}
		if (v1l < v2l)
			return v1l;
		else return v2l;
	}
}
