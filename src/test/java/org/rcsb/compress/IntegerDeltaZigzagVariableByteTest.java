package org.rcsb.compress;

import static org.junit.Assert.assertArrayEquals;

import java.util.Random;

import org.junit.Test;

public class IntegerDeltaZigzagVariableByteTest {

	@Test
	public void test1() {
		System.out.println("test1");
		int[] expecteds = {1,1,1,1,1,1,1,1,1,1,1,1,1};
		IntegerDeltaZigzagVariableByte transform = new IntegerDeltaZigzagVariableByte();
		byte[] out = transform.forward(expecteds);
		printInfo(expecteds, out);
		int[] actuals = transform.reverse(out);
		assertArrayEquals(expecteds, actuals);
	}

	@Test
	public void test2() {
		System.out.println("test2");
		int[] expecteds = {1,2,3,4,5,6,7,8,9,10,11,12,13};
		IntegerDeltaZigzagVariableByte transform = new IntegerDeltaZigzagVariableByte();
		byte[] out = transform.forward(expecteds);
		printInfo(expecteds, out);
		int[] actuals = transform.reverse(out);
		assertArrayEquals(expecteds, actuals);
	}
	
	@Test
	public void test3() {
		System.out.println("test3");
		int[] expecteds = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096};
		IntegerDeltaZigzagVariableByte transform = new IntegerDeltaZigzagVariableByte();
		byte[] out = transform.forward(expecteds);
		printInfo(expecteds, out);
		int[] actuals = transform.reverse(out);
		assertArrayEquals(expecteds, actuals);
	}
	
	@Test
	public void test4() {
		System.out.println("test4");
		int[] expecteds = {4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096};
		IntegerDeltaZigzagVariableByte transform = new IntegerDeltaZigzagVariableByte();
		byte[] out = transform.forward(expecteds);
		printInfo(expecteds, out);
		int[] actuals = transform.reverse(out);
		assertArrayEquals(expecteds, actuals);
	}
	
	@Test
	public void testRandom8Bit() {
		System.out.println("testRandom8Bit");
		Random r = new Random();
		int[] expecteds = new int[1000000];
		for (int i = 0; i < expecteds.length; i++) {
			expecteds[i] = r.nextInt(255) - 128;
		}
		
		long start = System.nanoTime();
		IntegerDeltaZigzagVariableByte transform = new IntegerDeltaZigzagVariableByte();
		byte[] out = transform.forward(expecteds);
		long end = System.nanoTime();
		long rate = (long) (expecteds.length/((end-start)*1E-9));
		System.out.println("Compression rate: [n=" + expecteds.length + "]: " + rate + " integers/second");
		
		printInfo(expecteds, out);
		
		start = System.nanoTime();
		int[] actuals = transform.reverse(out);
		end = System.nanoTime();
		rate = (long) (actuals.length/((end-start)*1E-9));
		System.out.println("Decompression rate: [n=" + expecteds.length + "]: " + rate + " integers/second");
		
		assertArrayEquals(expecteds, actuals);
	}
	
	@Test
	public void testRandom16() {
		System.out.println("testRandom16");
		Random r = new Random();
		int[] expecteds = new int[100];
		for (int i = 0; i < expecteds.length; i++) {
			expecteds[i] = r.nextInt(65536) - 32768;
		}
		
		IntegerDeltaZigzagVariableByte transform = new IntegerDeltaZigzagVariableByte();
		byte[] out = transform.forward(expecteds);
		printInfo(expecteds, out);
		int[] actuals = transform.reverse(out);
		assertArrayEquals(expecteds, actuals);
	}
	
	@Test
	public void testRandom32() {
		System.out.println("testRandom32");
		Random r = new Random();
		int[] expecteds = new int[100];
		for (int i = 0; i < expecteds.length; i++) {
			expecteds[i] = r.nextInt();
		}
		
		IntegerDeltaZigzagVariableByte transform = new IntegerDeltaZigzagVariableByte();
		byte[] out = transform.forward(expecteds);
		printInfo(expecteds, out);
		int[] actuals = transform.reverse(out);
		assertArrayEquals(expecteds, actuals);
	}
	
	private static void printInfo(int[] in, byte[] out) {
		System.out.println("Compress from: " + (in.length*4) + " to " + out.length + " bytes. Compression ratio: " + (in.length*4.0/out.length));
	}
}
