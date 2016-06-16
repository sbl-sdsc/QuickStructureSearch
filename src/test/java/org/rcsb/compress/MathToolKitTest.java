package org.rcsb.compress;

import java.util.Arrays;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import scala.util.Random;

public class MathToolKitTest {

	@Test
	public void test() {
		Random r = new Random();
	//	int expected = r.nextInt();
		int expected = 4096;
		System.out.println(Arrays.toString(MathToolKit.decompose(expected)));
		int[] multipliers = MathToolKit.decompose(expected);
		int actual = 0;
		for (int m: multipliers) {
			actual += Math.pow(2, m);
		}
	    assertEquals(expected, actual);
	}
}
