package org.rcsb.structuralSimilarity;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

public class SmithWatermanForIntHL implements PairFunction<Tuple2<Integer,Integer>,String,Float> {

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
    // index of max score
    private int maxX = 0;
    private int maxY = 0;
    // print traceback if it is greater than 0
    private int traceback = 0;


	public SmithWatermanForIntHL(Broadcast<List<Tuple2<String,Vector>>> data, int traceback) {
		this.data = data;
		this.traceback = traceback;
	}
	
	@Override
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) {
		Tuple2<String,Vector> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Vector> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		int[] v1 = toIntVector(t1._2.toArray());
		int[] v2 = toIntVector(t2._2.toArray());
		
		v1length = v1.length;
		v2length = v2.length;
		score = new double[v1length+1][v2length+1];
		score[0][0] = 0;
		int[][] b = new int[v1.length][v2.length];
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
		    	double maxScore = maximum(diagScore, upScore,leftScore, 0);
		    	if (maxScore == 0) {
		    		b[i-1][j-1] = -1;
		    	} else if (maxScore == diagScore) {
		    		b[i-1][j-1] = 0;
		    	} else if (maxScore == upScore) {
		    		b[i-1][j-1] = 1;
		    	} else {
		    		b[i-1][j-1] = 2;
		    	}
		    	score[i][j] = maxScore;
		    }
		}
		float value = (float) getMaxScore();
		if (v1length > v2length)
			value = value/v2length;
		else value = value/v1length;
		if (traceback > 0)  {
			printTraceback(v1,v2,b);
			System.out.println("value: "+value + " v1length: "+ v1length + " v2length: " + v2length);
			System.out.println();
		}
		return new Tuple2<String, Float>(key.toString(), value);
	}
	
	private void printTraceback(int[] v1,int[] v2,int[][] b) {
		int x = maxX;
		int y = maxY;
		String commonAngleV1 = " end at " + x;
		String commonAngleV2 = " end at " + y;
		while (x >= 0 && y >= 0 && b[x][y] >= 0) {
			if (b[x][y] == 0) {
				commonAngleV1 = String.format("% 3d",v1[x])+" "+commonAngleV1;
				commonAngleV2 = String.format("% 3d",v2[y])+" "+commonAngleV2;
				x--;
				y--;
			} else if (b[x][y] == 1) {
				commonAngleV1 = "--- "+commonAngleV1;
				commonAngleV2 = "    "+commonAngleV2;
				x--;
			} else if (b[x][y] == 2) {
				commonAngleV1 = "    "+commonAngleV1;
				commonAngleV2 = "--- "+commonAngleV2;
				y--;
			} else if (b[x][y] == -1) {
				break;
			}
		}
		x++;
		y++;
		commonAngleV1 = "start from "+ x + "\t" +commonAngleV1;
		commonAngleV2 = "start from "+ y + "\t" +commonAngleV2;
		System.out.println(commonAngleV1);
		System.out.println(commonAngleV2);
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
		maxX = 0;
		maxY = 0;
		// skip the first row and column
		for (int i = 1; i <= v1length; i++) {
		    for (int j = 1; j <= v2length; j++) {
				if (score[i][j] > maxScore) {
					maxX = i;
					maxY = j;
				    maxScore = score[i][j];
				}
		    }
		}
		maxX--;
		maxY--;
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

	private int[] toIntVector(double[] vector) {
		int[] intVector = new int[vector.length];
		for (int i = 0; i < vector.length; i++) {
			intVector[i] = (int)Math.round(vector[i]);
		}
		return intVector;
	}
}
