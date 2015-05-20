package org.rcsb.structuralSimilarity;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the length longest common subsequence over the length of the chains
 */
public class LCSFeatureIndexHL implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	// the amount to consider two value different
	private double diff = 3.14/72;
	private Broadcast<List<Tuple2<String,Vector>>> data = null;
	// print traceback if it is greater than 0
	private int traceback = 0;
	/**
     * Constructor with default setting
     */
	public LCSFeatureIndexHL(Broadcast<List<Tuple2<String,Vector>>> data,int traceback) {
		this.data = data;
		this.traceback = traceback;
	}
	
	/**
     * Constructor with parameter
     */
	public LCSFeatureIndexHL(Broadcast<List<Tuple2<String,Vector>>> data, double diff) {
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
	private int LCS(double[] v1, double[] v2, double diff) {
		// the length of LCS at each pair of values
		int c[][] = new int[v1.length][v2.length];
		// traceback
		int b[][] = new int[v1.length][v2.length];
		// the first row
		for (int i = 0; i < v1.length; i++) {
			if (similar(v1[i],v2[0],diff)) {
				c[i][0] = 1;
				b[i][0] = 0;
			} else {
				c[i][0] = 0;
				b[i][0] = -1;
			}
		}
		// the first column
		for (int j = 0; j < v2.length; j++) {
			if (similar(v1[0],v2[j],diff)) {
				c[0][j] = 1;
				b[0][j] = 0;
			}
			else {
				c[0][j] = 0;
				b[0][j] = -1;
			}
		}
		// the rest
		for (int i = 1; i < v1.length; i++) {
			for (int j = 1; j < v2.length; j++) {
				if (similar(v1[i],v2[j],diff)) {
					c[i][j] = c[i-1][j-1] + 1;
					b[i][j] = 0;
				}
				else if (c[i-1][j] >= c[i][j-1]) {
					c[i][j] = c[i-1][j];
					b[i][j] = 1;
				}
				else {
					c[i][j] = c[i][j-1];
					b[i][j] = 2;
				}
			}
		}
		if (traceback > 0)
			printLCS(v1,v2,b);
		// the last pair represent the global LCS
		return c[v1.length-1][v2.length-1];
	}
	
	private void printLCS(double[] v1, double[] v2, int[][] b) {
		int x = v1.length - 1;
		int y = v2.length - 1;
		//double[] commonAngleV1 = new double[v1.length];
		//double[] commonAngleV2 = new double[v2.length];
		while (x >= 0 && y >= 0) {
			if (b[x][y] == 0) {
				break;
			} else if (b[x][y] == 1) {
				x--;
			} else if (b[x][y] == 2) {
				y--;
			} else if (b[x][y] == -1) {
				break;
			}
		}
		String commonAngleV1 = " end at " + x;
		String commonAngleV2 = " end at " + y;
		while (x >= 0 && y >= 0) {
			if (b[x][y] == 0) {
				//commonAngleV1[x] = v1[x];
				//commonAngleV2[y] = v2[y];
				commonAngleV1 = String.format("%.2f",v1[x])+" "+commonAngleV1;
				commonAngleV2 = String.format("%.2f",v2[y])+" "+commonAngleV2;
				x--;
				y--;
			} else if (b[x][y] == 1) {
				commonAngleV1 = "---- "+commonAngleV1;
				commonAngleV2 = "     "+commonAngleV2;
				x--;
			} else if (b[x][y] == 2) {
				commonAngleV1 = "     "+commonAngleV1;
				commonAngleV2 = "---- "+commonAngleV2;
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
		System.out.println();
	}
	
	/**
	 * Compare the value
	 */
	private boolean similar(double x, double y, double diff) {
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
	private int getLength(double[] v1, double[] v2) {
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
