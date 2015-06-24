package org.rcsb.structuralSimilarity;

import java.util.List;

import org.rcsb.utils.BlastClustReader;

public class Demo_pr {

	public static void main(String[] args) {
		
		int sequenceIdentity = 100;
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);
		
		List<List<String>> clusters = reader.getPdbChainIdClusters();
		
		for (List<String> cluster: clusters) {
			if (cluster.size() == 6) {
				System.out.println(cluster);
			}
		}
	}
}
