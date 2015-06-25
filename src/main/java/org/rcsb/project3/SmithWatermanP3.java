package org.rcsb.project3;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the longest local common subsequence over the length of the chains
 * using the SmithWaterman algorithm
 * 
 * @author Chris Li
 */
public class SmithWatermanP3 implements PairFunction<Tuple2<Integer,Integer>,String,Float> {

	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data = null;
    // score of indel
    private double indel = 1;
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

    public SmithWatermanP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
    
    /**
     * Constructor with traceback option
     * @param data
     * @param traceback
     */
	public SmithWatermanP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data, int traceback) {
		this.data = data;
		this.traceback = traceback;
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
		
		v1length = v1.length();
		v2length = v2.length();
		score = new double[v1length+1][v2length+1];
		score[0][0] = 0;
		int[][] b = new int[v1length][v2length];
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
		    	double diagScore = score[i - 1][j - 1] + calSimilarity(v1,v2,i-1,j-1);
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
		// Traceback
		if (traceback > 0)  {
			printTraceback(v1,v2,b);
			System.out.println("value: "+value + " v1length: "+ v1length + " v2length: " + v2length);
			System.out.println();
		}
		return new Tuple2<String, Float>(key.toString(), value);
	}
	
	/**
	 * Calculate the similarity. Object class casting.
	 * @param v1
	 * @param v2
	 * @param i
	 * @param j
	 * @return
	 */
	@SuppressWarnings("unchecked")
	private <T,K> double calSimilarity(SequenceFeatureInterface<T> v1,SequenceFeatureInterface<K> v2,int i,int j) {
		return v1.similarity((SequenceFeatureInterface<T>)v2,i,j);
	}
	
	/**
	 * Print the SmithWaterman traceback
	 * @param v1
	 * @param v2
	 * @param b
	 */
	private void printTraceback(SequenceFeatureInterface<?> v1,SequenceFeatureInterface<?> v2,int[][] b) {
		int x = maxX;
		int y = maxY;
		String commonV1 = " end at " + x;
		String commonV2 = " end at " + y;
		while (x >= 0 && y >= 0 && b[x][y] >= 0) {
			if (b[x][y] == 0) {
				commonV1 = v1.toString(x)+" \t"+commonV1;
				commonV2 = v2.toString(y)+" \t"+commonV2;
				x--;
				y--;
			} else if (b[x][y] == 1) {
				commonV1 = "xxx \t"+commonV1;
				commonV2 = "--- \t"+commonV2;
				x--;
			} else if (b[x][y] == 2) {
				commonV1 = "--- \t"+commonV1;
				commonV2 = "xxx \t"+commonV2;
				y--;
			} else if (b[x][y] == -1) {
				break;
			}
		}
		if (x < 0)
			x++;
		if (y < 0)
			y++;
		commonV1 = "start from "+ x + "\t" +commonV1;
		commonV2 = "start from "+ y + "\t" +commonV2;
		System.out.println(commonV1);
		System.out.println(commonV2);
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
	
    /**
     * Get the maximum value of the inputs
     * @param a
     * @param b
     * @param c
     * @param d
     * @return
     */
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
