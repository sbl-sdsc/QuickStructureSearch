package org.rcsb.structuralSimilarity;

import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

/**
 * Splits a comma separated set of items at the n-th comma
 * 
 * @author Peter Rose
 *
 */
public class SplitAtComma implements PairFunction<String, String, String> {
	private static final long serialVersionUID = 1L;
	private int n;
	
	public SplitAtComma(int n) {
		this.n = n;
	}

	@Override
	public Tuple2<String, String> call(String t) throws Exception {
		int split = 0;
		for (int i = 0, commas = 0; i < t.length(); i++) {
			if (t.charAt(i) == ',') {
				commas++;
				if (commas == n) {
					split = i;
					break;
				}
			}
		}

		String key = t.substring(0, split);
		String value = t.substring(split+1);
		return new Tuple2<String,String>(key, value);
	}
}
