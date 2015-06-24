package org.rcsb.project3;

public class LevenshteinDistanceP3 {
	// http://en.wikipedia.org/wiki/Levenshtein_distance
	
	public static <T> int distance(SequenceFeatureInterface<T> s, SequenceFeatureInterface<T> t)
	{
	    if (s.length() == 0) return t.length();
	    if (t.length() == 0) return s.length();
	 
	    // create two work vectors of integer distances
	    int[] v0 = new int[t.length() + 1];
	    int[] v1 = new int[t.length() + 1];
	 
	    // initialize v0 (the previous row of distances)
	    // this row is A[0][i]: edit distance for an empty s
	    // the distance is just the number of characters to delete from t
	    for (int i = 0; i < v0.length; i++)
	        v0[i] = i;
	 
	    for (int i = 0; i < s.length(); i++) {
	        // calculate v1 (current row distances) from the previous row v0
	 
	        // first element of v1 is A[i+1][0]
	        //   edit distance is delete (i+1) chars from s to match empty t
	        v1[0] = i + 1;
	 
	        // use formula to fill in the rest of the row
	        for (int j = 0; j < t.length(); j++) {
	            int cost = s.identity(t, i, j) ? 0 : 1;
	            v1[j + 1] = Math.min(Math.min(v1[j] + 1, v0[j + 1] + 1), v0[j] + cost);
		    }
	 
	        // copy v1 (current row) to v0 (previous row) for next iteration
	        for (int j = 0; j < v0.length; j++)
	            v0[j] = v1[j];
	    }
	 
	    return v1[t.length()];
	}
	
	/**
	 * Returns value between 0 and 1
	 * @param s
	 * @param t
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public static <T,K> double normalizedDistance(SequenceFeatureInterface<T> s1, 
			SequenceFeatureInterface<K> s2) {
		int editDistance = distance(s1, (SequenceFeatureInterface<T>)s2);
		return 1 - editDistance/(double)Math.max(s1.length(),  s2.length());
	}
	
}
