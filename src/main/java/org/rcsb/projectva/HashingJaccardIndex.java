package org.rcsb.projectva;

import java.util.Map;
import java.util.Map.Entry;

/**
 * Class that calculates the Jaccard index between two feature maps.
 * @author Peter Rose
 *
 */
public class HashingJaccardIndex {

	/**
	 * Returns the Jaccard index for two feature maps. The feature maps represent multi-set,
	 * where elements can occur more than once. The Jaccard index is a number between 0 and 1.
	 * If the two feature maps are identical, the index is 1, if they share some features, the
	 * index will be between 0 and 1, and if no features are in common, the index will be 0.
	 * Each feature map consists of a key, that represents the feature of type T, 
	 * and a value that represent the feature count.
	 * @param featureCounts1 feature map
	 * @param featureCounts2 feature map
	 * @return the Jaccard index
	 */
	


	
	public static<T> double distance(int seq1[], int seq2[]) 
	{
		//seq1 and seq2 are necessarily the same size
		int min = 0;
		int max = 0;
		for (int i = 0; i < seq1.length; i++) {
			min += (int)Math.min(seq1[i], seq2[i]);
			max += (int)Math.max(seq1[i], seq2[i]);
		}
		
		double value = 0.0;
		if (max != 0) {
			value = (double)min/(double)max;
		}
		
		return value;
	}
		
		
}
