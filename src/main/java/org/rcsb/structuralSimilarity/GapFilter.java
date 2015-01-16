package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.Function;

import scala.Tuple2;

public class GapFilter implements Function<Tuple2<String,Point3d[]>, Boolean> {
	private static final long serialVersionUID = 1L;
	int maxGapSize;
	int maxGapCount;

	public GapFilter(int maxGapSize, int maxGapCount) {
		this.maxGapSize = maxGapSize;
		this.maxGapCount = maxGapCount;
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
		
		// scan the truncated chain for gaps
		int gapSize = 0;
		int gapCount = 0;
		boolean hasGap = false;
		
		for (int i = start; i < end; i++) {
			if (points[i] == null) {
				gapSize++;
				hasGap = true;
			} else {
				gapSize = 0;
				if (hasGap) {
					gapCount++;
					hasGap = false;
				}
			}
			if (gapSize > this.maxGapSize || gapCount > this.maxGapCount) {
				return false;
			}
		}
		return true;
	}
}
