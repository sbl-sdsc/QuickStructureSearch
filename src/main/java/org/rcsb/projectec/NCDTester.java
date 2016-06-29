package org.rcsb.projectec;

import org.rcsb.structuralSimilarity.NormalizedCompressionDistance;

public class NCDTester {
	public static void main(String[] args) {
		int[] a = {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5};
		int[] b = {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5};
		System.out.println(NormalizedCompressionDistance.distance(a, b));
		
		int[] c = {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5};
		int[] d = {5,4,3,2,1,5,4,3,2,1,5,4,3,2,1,5,4,3,2,1,5,4,3,2,1};
		System.out.println(NormalizedCompressionDistance.distance(c, d));
		
		int[] e = {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5};
		int[] f = {1,3,4,4,5,1,2,3,4,5,1,5,3,4,5,1,2,7,4,5,1,2,3,4,5};
		System.out.println(NormalizedCompressionDistance.distance(e, f));
		
		int[] g = {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5};
		int[] h = {7,8,9,6,7,8,5,4,3,6,7,8,9,7,5,7,5,7,2,5,9,9,9,9,9,9,9,9,78,54,35,43543,54,3,43243,423,4,32,432432432,4};
		System.out.println(NormalizedCompressionDistance.distance(g, h));
	}
}
