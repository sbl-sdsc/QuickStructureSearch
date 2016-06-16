package org.rcsb.compress;

import java.util.Random;

public class TestDataSets {

	private static int LENGTH = 3000;
	private static int SEED = 1;
	
	private static int[] random8Bit = new int[LENGTH];
	private static int[] random16Bit = new int[LENGTH];
	private static int[] random32Bit = new int[LENGTH];
	
	static {
		Random r = new Random(SEED);

		for (int i = 0; i < LENGTH; i++) {
			random8Bit[i] = r.nextInt(255) - 128;
			random16Bit[i] = r.nextInt(65536) - 32768;
			random32Bit[i] = r.nextInt();
		}
	}
	
	private static int[][] testCases = {
//			{1,1,1,1,1,1,1,1,1,1,1,1},
			{1,2,3,4,5,6,7,8,9,10,11,12},
			{1,2,4,8,16,32,64,128,256,512,1024,2048},
			{4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096},
			random8Bit,
			random16Bit,
			random32Bit
	};
	
	public static int [][] getTestCases() {
		return testCases;
	}
}
