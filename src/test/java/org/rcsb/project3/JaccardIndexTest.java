package org.rcsb.project3;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

public class JaccardIndexTest {

	@Test
	public void testIdentity() {
		Map<Integer, Integer> features1 = new HashMap<Integer, Integer>();
		features1.put(1,1);
		features1.put(2,1);
		features1.put(3,1);
		
		Map<Integer, Integer> features2 = new HashMap<Integer, Integer>();
		features2.put(1,1);
		features2.put(2,1);
		features2.put(3,1);
		
		double index = JaccardIndex.jaccardIndex(features1, features2);
		assertEquals(1.0, index, Math.ulp(1.0));
	}
	
	@Test
	public void testIdentity2() {
		Map<Integer, Integer> features1 = new HashMap<Integer, Integer>();
		features1.put(1,2);
		features1.put(2,2);
		features1.put(3,2);
		
		Map<Integer, Integer> features2 = new HashMap<Integer, Integer>();
		features2.put(1,2);
		features2.put(2,2);
		features2.put(3,2);
		
		double index = JaccardIndex.jaccardIndex(features1, features2);
		assertEquals(1.0, index, Math.ulp(1.0));
	}
	
	@Test
	public void testNonIdentity() {
		Map<Integer, Integer> features1 = new HashMap<Integer, Integer>();
		features1.put(1,1);
		features1.put(2,1);
		features1.put(3,1);
		
		Map<Integer, Integer> features2 = new HashMap<Integer, Integer>();
		features2.put(4,1);
		features2.put(5,1);
		features2.put(6,1);
		
		double index = JaccardIndex.jaccardIndex(features1, features2);
		assertEquals(0.0, index, Math.ulp(1.0));
	}
	
	@Test
	public void testUnequalFeatureCounts() {
		Map<Integer, Integer> features1 = new HashMap<Integer, Integer>();
		features1.put(1,1);
		features1.put(2,1);
		features1.put(3,1);
		
		Map<Integer, Integer> features2 = new HashMap<Integer, Integer>();
		features2.put(1,2);
		features2.put(2,2);
		features2.put(3,2);
		
		double index = JaccardIndex.jaccardIndex(features1, features2);
		assertEquals(0.5, index, Math.ulp(1.0));
	}
	
	@Test
	public void testUnequalFeatureCounts2() {
		Map<Integer, Integer> features1 = new HashMap<Integer, Integer>();
		features1.put(1,1);
		features1.put(2,1);
		features1.put(3,1);
		
		Map<Integer, Integer> features2 = new HashMap<Integer, Integer>();
		features2.put(1,4);
		features2.put(2,4);
		features2.put(3,4);
		
		double index = JaccardIndex.jaccardIndex(features1, features2);
		assertEquals(0.25, index, Math.ulp(1.0));
	}
}
