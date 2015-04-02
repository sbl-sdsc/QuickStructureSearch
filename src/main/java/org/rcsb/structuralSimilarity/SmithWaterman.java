package org.rcsb.structuralSimilarity;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

public class SmithWaterman implements PairFunction<Tuple2<Integer,Integer>,String,Float> {

	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,Vector>>> data = null;
	// the amount to consider two value different
	private double diff = 3.14/72;
	// score of match
	private double match = 1;
	// score of mismatch
    private double mismatch = -1;
    // score of indel
    private double indel = -1;
    // score of gap
    private double gap = -1;
    // length of v1
    private int v1length;
    // length of v2
    private int v2length;
    // score matrix
    private double [][] score;


	public SmithWaterman(Broadcast<List<Tuple2<String,Vector>>> data) {
		this.data = data;
	}
	
	@Override
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) {
		Tuple2<String,Vector> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Vector> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		double[] v1 = t1._2.toArray();
		double[] v2 = t2._2.toArray();
		
		v1length = v1.length;
		v2length = v2.length;
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
		for (int i = 1; i <= v1length; i++) {
		    for (int j = 1; j <= v2length; j++) {
		    	double diagScore = score[i - 1][j - 1] + similarity(v1[i-1], v2[j-1]);
		    	double upScore = score[i][j - 1] + indel;
		    	double leftScore = score[i - 1][j] + indel;
		    	score[i][j] = maximum(diagScore, upScore,leftScore, 0);
		    }
		}
		float value = (float) getMaxScore();
		if (v1length > v2length)
			value = value/v2length;
		else value = value/v1length;
		return new Tuple2<String, Float>(key.toString(), value);
	}
	
	private double similarity(double x, double y) {
		if (Double.isNaN(x) || Double.isNaN(y)){
			return gap;
		}
		if (((x+diff) > y) && ((x-diff) < y))
			return match;
		else return mismatch;
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
	
	private static double maximum(double a, double b, double c, double d) {
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

}
