package org.rcsb.projectm;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

public class NCDIndexTest {

	@Test
	public void testIdentity() {
		int[] features1 = new int[25];
		for(int i = 1; i < 26;i++) 
		{
			features1[i - 1] = i;
		}
		
		int[] features2 = new int[25];
		for(int i = 1; i < 26;i++) 
		{
			features2[i - 1] = i;
		}
		
		double index = NormalizedCompressionDistance.distance(features1, features2);
		System.out.println(index);
		//assertEquals(0, index, Math.ulp(1.0));
	}
	
	@Test
	public void testReverse() {
		int[] features1 = new int[25];
		for(int i = 1; i < 26;i++) 
		{
			features1[i - 1] = i;
		}
		
		int[] features2 = new int[25];
		for(int i = 25; i >= 1;i--) 
		{
			features2[25 - i] = i;
		}
		
		double index = NormalizedCompressionDistance.distance(features1, features2);
		System.out.println(index);
		//assertEquals(0, index, Math.ulp(1.0));
	}
	
	@Test
	public void testSimilar() {
		int[] features1 = new int[25];
		for(int i = 1; i < 26;i++) 
		{
			features1[i - 1] = i;
		}
		
		int[] features2 = new int[25];
		for(int i = 1; i < 26;i++) 
		{
			features2[i - 1] = i + 4;
		}
		
		double index = NormalizedCompressionDistance.distance(features1, features2);
		System.out.println(index);
		//assertEquals(0, index, Math.ulp(1.0));
	}
	
	@Test
	public void testDissimilar() {
		int[] features1 = new int[25];
		for(int i = 1; i < 26;i++) 
		{
			features1[i - 1] = i;
		}
		
		int[] features2 = new int[25];
		for(int i = 1; i < 26;i++) 
		{
			features2[i - 1] = i + 100;
		}
		
		double index = NormalizedCompressionDistance.distance(features1, features2);
		System.out.println(index);
		//assertEquals(0, index, Math.ulp(1.0));
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
		assertEquals(1.0, index, Math.ulp(1.0));
	}
}
