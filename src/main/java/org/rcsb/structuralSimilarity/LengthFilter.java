package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.Function;

import scala.Tuple2;

/**
 * This class implements a filter in Spark workflow to filter out protein chains
 * that are either short than the specified minimum length or longer than the
 * specified maximum length. It ignores N-terminal and C-terminal gaps.
 * 
 * @author Peter Rose
 */
public class LengthFilter implements Function<Tuple2<String,Point3d[]>, Boolean> {
	private static final long serialVersionUID = 1L;
	int minLength;
	int maxLength;

	/**
	 * @param minLength minimum chain length
	 * @param maxLength maximum chain length
	 */
	public LengthFilter(int minLength, int maxLength) {
		this.minLength = minLength;
		this.maxLength = maxLength;
	}

	@Override
	public Boolean call(Tuple2<String,Point3d[]> tuple) throws Exception {
		Point3d[] points = tuple._2;
		
		int start = 0;
		
		// skip N-terminal gap (start of chain)
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				start = i;
				break;
			} 
		}
		
		// skip C-terminal gap
		int end = points.length-1;
		for (int i = points.length-1; i > start; i--) {
			if (points[i] != null) {
				end = i;
				break;
			}
		}
		
		int length = end - start + 1;
		
		return (length >= this.minLength && length <= this.maxLength);
	}
}
