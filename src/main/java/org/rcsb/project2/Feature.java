package org.rcsb.project2;

/**
 * Interface for a protein feature.
 * 
 * @author Kevin Wu
 *
 */
public interface Feature {
	/**
	 * Tests whether the ith set of distances match a given feature. This is used to find where a feature starts
	 * 
	 * @param i
	 *            integer for the distance between how many points
	 * @param d
	 *            the distance
	 * @return true if it matches the feature
	 */
	boolean match(int i, double d);

	/**
	 * For continuing matches
	 * 
	 * @param d
	 *            Array of double for distances
	 * @return true if it matches the feature
	 */
	boolean match(double[] d);
}
