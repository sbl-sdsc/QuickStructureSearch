package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.Function;

import scala.Tuple2;

public class LengthFilter implements Function<Tuple2<String,Point3d[]>, Boolean> {
	private static final long serialVersionUID = 1L;
	int minLength;
	int maxLength;

	public LengthFilter(int minLength, int maxLength) {
		this.minLength = minLength;
		this.maxLength = maxLength;
	}

	@Override
	public Boolean call(Tuple2<String,Point3d[]> tuple) throws Exception {
		Point3d[] points = tuple._2;
		
		int start = 0;
		
		// skip N-terminal gap
		for (int i = 0; i < points.length; i++) {
			if (points[i] == null) {
				start = i;
			} else {
				break;
			}
		}
		
		// skip C-terminal gap
		int end = points.length-1;
		for (int i = points.length-1; i > start; i--) {
			if (points[i] == null) {
				end = i;
			} else {
				break;
			}
		}
		
		int length = end - start + 1;
		
		return (length >= this.minLength && length <= this.maxLength);
	}
}
