package org.rcsb.structuralSimilarity;

import java.util.Arrays;


public class LevenshteinDistance {
	// http://en.wikipedia.org/wiki/Levenshtein_distance

	
	public static void main(String[] args) {
		int[] v1 = {1,2,3,4,5,6,7,8,9};
		int[] v2 = {2,3,4,6};
		System.out.println(distance(v1, v2));
		System.out.println(normalizedDistance(v1, v2));
		System.out.println(distanceOpt(v1, v2, 0.4));
		System.out.println(damerauLevenshteinDistance(v1, v2, 0.439));
		String[] alignment = alignG("ABCDEF","BCDXXF");
		System.out.println(alignment[0]);
		System.out.println(alignment[1]);
		alignment = align(v1,v2);
		System.out.println(alignment[0]);
		System.out.println(alignment[1]);
	}
	
	public static int distance(int[] s, int[] t)
	{
	    if (s.length == 0) return t.length;
	    if (t.length == 0) return s.length;
	 
	    // create two work vectors of integer distances
	    int[] v0 = new int[t.length + 1];
	    int[] v1 = new int[t.length + 1];
	 
	    // initialize v0 (the previous row of distances)
	    // this row is A[0][i]: edit distance for an empty s
	    // the distance is just the number of characters to delete from t
	    for (int i = 0; i < v0.length; i++)
	        v0[i] = i;
	 
	    for (int i = 0; i < s.length; i++) {
	        // calculate v1 (current row distances) from the previous row v0
	 
	        // first element of v1 is A[i+1][0]
	        //   edit distance is delete (i+1) chars from s to match empty t
	        v1[0] = i + 1;
	 
	        // use formula to fill in the rest of the row
	        for (int j = 0; j < t.length; j++) {
	            int cost = s[i] == t[j] ? 0 : 1;
	            v1[j + 1] = Math.min(Math.min(v1[j] + 1, v0[j + 1] + 1), v0[j] + cost);
		    }
	 
	        // copy v1 (current row) to v0 (previous row) for next iteration
	        for (int j = 0; j < v0.length; j++)
	            v0[j] = v1[j];
//	        System.out.println(v1[t.length]);
	    }
	 
//	    System.out.println("v1: " + Arrays.toString(v1));
	    return v1[t.length];
	}
	
	public static double distanceOpt(int[] source, int[] target, double threshold)
	{
		int s[] = source;
		int t[] = target;
		if (source.length > target.length) {
			s = target;
			t = source;
		}
	    if (s.length == 0) return 0;
	    if (t.length == 0) return 0;
	    
	    if ((1 - Math.abs(s.length - t.length)/(double) t.length) < threshold) { return 0.01; }
	 
	    // create two work vectors of integer distances
	    int[] v0 = new int[t.length + 1];
	    int[] v1 = new int[t.length + 1];
	 
	    // initialize v0 (the previous row of distances)
	    // this row is A[0][i]: edit distance for an empty s
	    // the distance is just the number of characters to delete from t
	    for (int i = 0; i < v0.length; i++)
	        v0[i] = i;
	 
	    for (int i = 0; i < s.length; i++) {
	        // calculate v1 (current row distances) from the previous row v0
	 
	        // first element of v1 is A[i+1][0]
	        //   edit distance is delete (i+1) chars from s to match empty t
	        v1[0] = i + 1;
	 
	        // use formula to fill in the rest of the row
	        for (int j = 0; j < t.length; j++) {	        
	           int cost = s[i] == t[j] ? 0 : 1;
	           v1[j + 1] = Math.min(Math.min(v1[j] + 1, v0[j + 1] + 1), v0[j] + cost);
		    }
	 
	        // early termination (doesn't work!!)
	//        System.out.println("minDist:: " + v1[t.length]);
	////        if ((1.0 - v1[t.length]/(double) t.length) < threshold) return  0;
	        
	        // copy v1 (current row) to v0 (previous row) for next iteration
	        for (int j = 0; j < v0.length; j++)
	            v0[j] = v1[j];
//	        System.out.println(v1[t.length]);
	    }
	 
//	    System.out.println("v1: " + Arrays.toString(v1));
	//    return v1[t.length];
	    return 1 - v1[t.length]/(double) t.length;
	}
	
	
	public static String[] align(String a, String b) {
        int[][] T = new int[a.length() + 1][b.length() + 1];

        for (int i = 0; i <= a.length(); i++)
            T[i][0] = i;

        for (int i = 0; i <= b.length(); i++)
            T[0][i] = i;

        for (int i = 1; i <= a.length(); i++) {
            for (int j = 1; j <= b.length(); j++) {
                if (a.charAt(i - 1) == b.charAt(j - 1))
                    T[i][j] = T[i - 1][j - 1];
                else
                    T[i][j] = Math.min(T[i - 1][j], T[i][j - 1]) + 1;
            }

        }

        StringBuilder aa = new StringBuilder(), bb = new StringBuilder();

        for (int i = a.length(), j = b.length(); i > 0 || j > 0; ) {
            if (i > 0 && T[i][j] == T[i - 1][j] + 1) {
                aa.append(a.charAt(--i));
                bb.append("-");
            } else if (j > 0 && T[i][j] == T[i][j - 1] + 1) {
                bb.append(b.charAt(--j));
                aa.append("-");
            } else if (i > 0 && j > 0 && T[i][j] == T[i - 1][j - 1]) {
                aa.append(a.charAt(--i));
                bb.append(b.charAt(--j));
            }
        }

        return new String[]{aa.reverse().toString(), bb.reverse().toString()};
    }
	
	public static String[] alignG(String a, String b) {
        int[][] T = new int[a.length() + 1][b.length() + 1];

        for (int i = 0; i <= a.length(); i++)
            T[i][0] = i;

        for (int i = 0; i <= b.length(); i++)
            T[0][i] = i;

        for (int i = 1; i <= a.length(); i++) {
            for (int j = 1; j <= b.length(); j++) {
                if (a.charAt(i - 1) == b.charAt(j - 1))
                    T[i][j] = T[i - 1][j - 1];
                else
                    T[i][j] = Math.min(T[i - 1][j], T[i][j - 1]) + 1;
            }

        }

        StringBuilder aa = new StringBuilder(), bb = new StringBuilder();
        
        boolean agap = false;
        boolean bgap = false;
        int gapPenalty = 0;

        for (int i = a.length(), j = b.length(); i > 0 || j > 0; ) {
            if (i > 0 && T[i][j] == T[i - 1][j] + 1) {
                aa.append(a.charAt(--i));
                bb.append("-");
                if (bgap) {
                	System.out.println("b gap extension");
                	gapPenalty+=1;
                } else {
                	System.out.println("new b gap");
                	gapPenalty+=10;
                	bgap = true;
                	agap = false;
                }
            } else if (j > 0 && T[i][j] == T[i][j - 1] + 1) {
                bb.append(b.charAt(--j));
                aa.append("-");
                if (agap) {
                	System.out.println("a gap extension");
                	gapPenalty+=1;
                } else {
                	System.out.println("new a gap");
                	gapPenalty+=10;
                	agap = true;
                	bgap = false;
                }
            } else if (i > 0 && j > 0 && T[i][j] == T[i - 1][j - 1]) {
                aa.append(a.charAt(--i));
                bb.append(b.charAt(--j));
                agap = false;
            }
        }
        System.out.println("gap penalty: " + gapPenalty);

        return new String[]{aa.reverse().toString(), bb.reverse().toString()};
    }
	
	public static String[] align(int[] a, int[] b) {
        int[][] T = new int[a.length + 1][b.length + 1];

        for (int i = 0; i <= a.length; i++)
            T[i][0] = i;

        for (int i = 0; i <= b.length; i++)
            T[0][i] = i;

        for (int i = 1; i <= a.length; i++) {
            for (int j = 1; j <= b.length; j++) {
                if (a[i - 1] == b[j - 1])
                    T[i][j] = T[i - 1][j - 1];
                else
                    T[i][j] = Math.min(T[i - 1][j], T[i][j - 1]) + 1;
            }

        }
        
        for (int i = 0; i < a.length+1; i++) {
            System.out.println(Arrays.toString(T[i]));
        }

        StringBuilder aa = new StringBuilder(), bb = new StringBuilder();

        for (int i = a.length, j = b.length; i > 0 || j > 0; ) {
            if (i > 0 && T[i][j] == T[i - 1][j] + 1) {
            	System.out.println("gap: " + T[i][j] + " " + i + "," + j + "->" + (i-1) + "," + j + " " + T[i - 1][j]);
                aa.append(a[--i]);
                bb.append("-");
            } else if (j > 0 && T[i][j] == T[i][j - 1] + 1) {
            	System.out.println("ins: " + T[i][j]  + " " + i + "," + j + "->" + (i) + "," + (j-1) + " " + T[i][j - 1]);
                bb.append(b[--j]);
                aa.append("-");
            } else if (i > 0 && j > 0 && T[i][j] == T[i - 1][j - 1]) {
            	System.out.println("mat: " + T[i][j]  + " " + i + "," + j + "->" + (i-1) + "," + (j-1) + " " + T[i - 1][j - 1]);
                aa.append(a[--i]);
                bb.append(b[--j]);
            }
        }

        return new String[]{aa.reverse().toString(), bb.reverse().toString()};
    }
	
	public static boolean lengthCheck(int[] s, int[] t, double threshold) {
		float maxLen = Math.max(s.length, t.length);
		float minLen = Math.min(s.length, t.length);
		return minLen/maxLen > threshold;
	}
	
	/**
	 * Returns value between 0 and 1
	 * @param s
	 * @param t
	 * @return
	 */
	public static double normalizedDistance(int[] s, int[] t) {
		int editDistance = distance(s, t);
		return 1 - editDistance/(double)Math.max(s.length,  t.length);
	}
	
	// http://stackoverflow.com/questions/9453731/how-to-calculate-distance-similarity-measure-of-given-2-strings
	/// Computes the Damerau-Levenshtein Distance between two strings, represented as arrays of
	/// integers, where each integer represents the code point of a character in the source string.
	/// Includes an optional threshhold which can be used to indicate the maximum allowable distance.
	/// </summary>
	/// <param name="source">An array of the code points of the first string</param>
	/// <param name="target">An array of the code points of the second string</param>
	/// <param name="threshold">Maximum allowable distance</param>
	/// <returns>Int.MaxValue if threshhold exceeded; otherwise the Damerau-Leveshteim distance between the strings</returns>
	public static double damerauLevenshteinDistance(int[] source, int[] target, double threshold) {
		int[] s = source;
		int[] t = target;
		
	    // Ensure arrays [i] / length1 use shorter length 
		if (source.length > target.length) {
			t = source;
			s = target;
		}

	    int length1 = s.length;
	    int length2 = t.length;

	    // Return trivial case - difference in string lengths exceeds threshhold
//	    if (Math.abs(length1 - length2) > threshold) { return Integer.MAX_VALUE; }
//	    if ((1 - Math.abs(length1 - length2)/(double) length2) < threshold) { return 0.01; }

	    // Ensure arrays [i] / length1 use shorter length 


	    int maxi = length1;
	    int maxj = length2;

	    int[] dCurrent = new int[maxi + 1];
	    int[] dMinus1 = new int[maxi + 1];
	    int[] dMinus2 = new int[maxi + 1];
	    int[] dSwap;

	    for (int i = 0; i <= maxi; i++) { dCurrent[i] = i; }

	    int jm1 = 0, im1 = 0, im2 = -1;

	    for (int j = 1; j <= maxj; j++) {

	        // Rotate
	        dSwap = dMinus2;
	        dMinus2 = dMinus1;
	        dMinus1 = dCurrent;
	        dCurrent = dSwap;

	        // Initialize
	        int minDistance = Integer.MAX_VALUE;
	        dCurrent[0] = j;
	        im1 = 0;
	        im2 = -1;

	        for (int i = 1; i <= maxi; i++) {

	            int cost = s[im1] == t[jm1] ? 0 : 1;

	            int del = dCurrent[im1] + 1;
	            int ins = dMinus1[i] + 1;
	            int sub = dMinus1[im1] + cost;

	            //Fastest execution for min value of 3 integers
	            int min = (del > ins) ? (ins > sub ? sub : ins) : (del > sub ? sub : del);

	            if (i > 1 && j > 1 && s[im2] == t[jm1] && s[im1] == t[j - 2])
	                min = Math.min(min, dMinus2[im2] + cost);

	            dCurrent[i] = min;
	            if (min < minDistance) { minDistance = min; }
	            im1++;
	            im2++;
	        }
	        jm1++;
//	        System.out.println("minDist: " + minDistance);
	//        if (minDistance > threshold) { return Integer.MAX_VALUE; }
	//        if ((1.0 - minDistance/(double) length2) < threshold) return  0;
	    }

//	    int result = dCurrent[maxi];
	    return 1 - dCurrent[maxi]/(double) length2;
//	    return (result > threshold) ? Integer.MAX_VALUE : result;
	}
	
}
