package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.Function;

import scala.Tuple2;
/**
* This class implements a filter in a Spark workflow to filter out protein chains
* that have at least one gap larger than the maximum specified gap size
* or that have more gaps than the maximum specified gap count.
* It ignores N-terminal and C-terminal gaps.
* 
* @author Peter Rose
*/
public class GapFilter implements Function<Tuple2<String,Point3d[]>, Boolean> {
	private static final long serialVersionUID = 1L;
	private int maxGapSize;
	private int maxGapCount;

	/**
	 * Constructor for GapFiler
	 * @param maxGapSize maximum gap size
	 * @param maxGapCount maximum gap count
	 */
	public GapFilter(int maxGapSize, int maxGapCount) {
		this.maxGapSize = maxGapSize;
		this.maxGapCount = maxGapCount;
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
		
		if (start == points.length-1) {
			return false;
		}
		
		// skip C-terminal gap (end of chain)
		int end = points.length-1;
		for (int i = points.length-1; i > start; i--) {
			if (points[i] != null) {
				end = i;
				break;
			}
		}
		
		// scan the truncated chain for gaps
		int gapSize = 0;
		int gapCount = 0;
		int maximumGapSize = 0;
		boolean hasGap = false;
		
		for (int i = start; i <= end; i++) {
			if (points[i] == null) {
				gapSize++;
				hasGap = true;
			} else {
				maximumGapSize = Math.max(gapSize, maximumGapSize);
				gapSize = 0;
				if (hasGap) {
					gapCount++;
					hasGap = false;
				}
			}
			maximumGapSize = Math.max(gapSize, maximumGapSize);
			if (maximumGapSize > this.maxGapSize || gapCount > this.maxGapCount) {
				return false;
			}
		}
		return true;
	}
}
