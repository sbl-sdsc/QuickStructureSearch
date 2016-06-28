package org.rcsb.rehsDavidM;

import java.util.Map;
import java.util.Map.Entry;

/**
 * Class that calculates the Jaccard index between two feature maps.
 * @author David Mao
 *
 */
public class MeetMinIndex {

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
	public static <T> double meetMinIndex(Map<T, Integer> featureCounts1, Map<T, Integer> featureCounts2) {
		int size1 = 0;
		int size2 = 0;
		int intersection = 0;
		
		for (Entry<T, Integer> entry: featureCounts1.entrySet()) {
			size1 += entry.getValue();
			if (featureCounts2.containsKey(entry.getKey())) {
				int v1 = entry.getValue();
				int v2 = featureCounts2.get(entry.getKey());
				intersection += Math.min(v1, v2);

			}
		}
		
		// Add remaining features that are exclusively in feature map 2
		for (Entry<T, Integer> entry: featureCounts2.entrySet()) {
			size2 += entry.getValue();
		}
		
		double value = 0;
		if (size1 > 0 && size2 > 0) {
		    value = intersection/(float)Math.min(size1, size2);
		}
		return value;
	}
}
