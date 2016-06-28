package org.rcsb.rehsDavidM;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;
/**
 * 
 * @author David Mao
 *
 */
public class NormalizedCompressionDistanceTest {

	@Test
	public void testIdentity() {
		int[] features1 = new int[26];
		for (int i = 0; i < 26; i++){
			features1[i] = i;
		}
		int[] features2 = new int[26];
		for (int i = 0; i < 26; i++){
			features2[i] = i;
		}
		double index = NormalizedCompressionDistance.distance(features1, features2);
		System.out.println(index);
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
		
		double index = MeetMinIndex.meetMinIndex(features1, features2);
		assertEquals(1.0, index, Math.ulp(1.0));
	}
	
	@Test
	public void testNonIdentity() {

		int[] features1 = new int[26];
		for (int i = 0; i < 26; i++){
			features1[i] = i;
		}
		int[] features2 = new int[26];
		for (int i = 0; i < 26; i++){
			features2[i] = 26-i;
		}
		double index = NormalizedCompressionDistance.distance(features1, features2);
		System.out.println(index);
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
		
		double index = MeetMinIndex.meetMinIndex(features1, features2);
		assertEquals(1, index, Math.ulp(1.0));
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
		
		double index = MeetMinIndex.meetMinIndex(features1, features2);
		assertEquals(1, index, Math.ulp(1.0));
	}
}
