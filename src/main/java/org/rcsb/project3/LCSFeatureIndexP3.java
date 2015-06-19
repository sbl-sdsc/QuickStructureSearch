package org.rcsb.project3;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the length longest common subsequence over the length of the chains
 * @author Chris Li
 */
public class LCSFeatureIndexP3 implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data = null;
	// print traceback if it is greater than 0
	private int traceback = 0;
	/**
     * Constructor with traceback option
     */
	public LCSFeatureIndexP3(Broadcast<List<Tuple2<String, SequenceFeatureInterface<?>>>> featureVectorsBc,int traceback) {
		this.data = featureVectorsBc;
		this.traceback = traceback;
	}
	
	public LCSFeatureIndexP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}

	/**
	 * Returns <PdbId.Chain, LCS score> pairs.
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,SequenceFeatureInterface<?>> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,SequenceFeatureInterface<?>> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		// get the 2 sequence
		SequenceFeatureInterface<?> v1 = t1._2;
		SequenceFeatureInterface<?> v2 = t2._2;
		
		// The longest common substring
		float LCSres = (float)LCS(v1,v2);
		float value = LCSres/(float)getLength(v1,v2);
        return new Tuple2<String, Float>(key.toString(), value);
    }
	
	/**
	 * Get the length of the longest common substring
	 * @param <K>
	 */
	@SuppressWarnings("unchecked")
	private <T, K> int LCS(SequenceFeatureInterface<T> v1, SequenceFeatureInterface<K> v2) {
		// the length of LCS at each pair of values
		int c[][] = new int[v1.length()][v2.length()];
		// traceback
		int b[][] = new int[v1.length()][v2.length()];
		// the first row
		for (int i = 0; i < v1.length(); i++) {
			if (v1.identity((SequenceFeatureInterface<T>)v2, i, 0)) {
				c[i][0] = 1;
				b[i][0] = 0;
			} else {
				c[i][0] = 0;
				b[i][0] = -1;
			}
		}
		// the first column
		for (int j = 0; j < v2.length(); j++) {
			if (v1.identity((SequenceFeatureInterface<T>)v2, 0, j)) {
				c[0][j] = 1;
				b[0][j] = 0;
			}
			else {
				c[0][j] = 0;
				b[0][j] = -1;
			}
		}
		// the rest
		for (int i = 1; i < v1.length(); i++) {
			for (int j = 1; j < v2.length(); j++) {
				if (v1.identity((SequenceFeatureInterface<T>)v2, i, j)) {
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
		// Traceback
		if (traceback > 0)
			printLCS(v1,v2,b);
		// the last pair represent the global LCS
		return c[v1.length()-1][v2.length()-1];
	}
	
	/**
	 * Print the LCS traceback
	 * @param v1
	 * @param v2
	 * @param b
	 */
	private void printLCS(SequenceFeatureInterface<?> v1, SequenceFeatureInterface<?> v2, int[][] b) {
		int x = v1.length() - 1;
		int y = v2.length() - 1;
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
		String commonV1 = " end at " + x;
		String commonV2 = " end at " + y;
		while (x >= 0 && y >= 0) {
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
		System.out.println();
	}
	
	/**
	 * Return the length of the shorter sequence
	 */
	private int getLength(SequenceFeatureInterface<?> v1, SequenceFeatureInterface<?> v2) {
		int v1l = v1.length();
		int v2l = v2.length();
		// gap is not counted
		for (int i = 0; i < v1.length(); i++) {
			if (v1.get(i) == null)
				v1l--;
		}
		for (int i = 0; i < v2.length(); i++) {
			if (v2.get(i) == null)
				v2l--;
		}
		if (v1l < v2l)
			return v1l;
		else return v2l;
	}
}
