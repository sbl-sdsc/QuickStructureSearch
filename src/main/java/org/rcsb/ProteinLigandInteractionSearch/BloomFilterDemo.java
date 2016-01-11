package org.rcsb.ProteinLigandInteractionSearch;

import java.util.Arrays;
import java.util.List;

import com.sangupta.bloomfilter.BloomFilter;
import com.sangupta.bloomfilter.impl.InMemoryBloomFilter;

public class BloomFilterDemo {

	public static void main(String[] args) {
	    int n = 100;
	    double falsePositiveRate = 0.01;
	    
	    List<String> sample = Arrays.asList("HisNE2AspOD2-24","HisND1AspOD2-24");
	    List<String> test = Arrays.asList("HisNE2AspOD2-24","HisND1AspOD1-24","HisCAAspOD1-24","HisND1AspOD2-23","HisND1AspOD2-24");
		
	    BloomFilter<String> filter = new InMemoryBloomFilter<String>(n, falsePositiveRate);
	    System.out.println("Bits: " + filter.getNumberOfBits());
        filter.addAll(sample);
		
        for (String example: test) {
        	System.out.println(example + ": " + filter.contains(example));
        }
	    
	}
}
