package org.rcsb.compress;

import static org.junit.Assert.assertArrayEquals;

import java.util.Arrays;

import org.junit.Test;
import org.rcsb.compress.dev.Ivqx5Transform;

public class TransformTest {

	private IntegerTransform[] transforms = { 
			new NoOpTransform(),
			new ColToRowTransform(),
			new CombinedTransform(new NoOpTransform(), new ColToRowTransform()),
			new Delta4Transform(),
			new DeltaTransform(),
			new DiagonalXorTransform(),
			new HorizontalXorTransform(),
			new AncientEgyptianDecomposition(new LeGallWavelet()),
			new RefDeltaTransform(),
			new UnsignedDeltaTransform(),
			new UnsignedRefDeltaTransform(),
			new UnsignedTransform(),
			new XorTransform()
//			new Ivqx5Transform()

//			new SphericalCoordinateTransform() // fails on test 3 ...
//			new SixSphereCoordinateTransform() // fails on test 5 ...
//			new AncientEgyptianDecomposition(new IntCompressorTransform(new DeltaZigzagVariableByte()))
//			new AncientEgyptianDecomposition(new FastWaveletTransform("Daubechies 2") // fails for all but test 0-1)
//			new AncientEgyptianDecomposition(new FastWaveletTransform("Haar orthogonal") // fails for random32Bit set)
			};
	
	private int[][] testSets = TestDataSets.getTestCases();

	@Test
	public void test() {
		for (IntegerTransform t : transforms) {

			for (int i = 0; i < testSets.length; i++) {
				System.out.println(t + ": testset: " + i);
				int[] expecteds = testSets[i];
				int[] out = t.forward(expecteds);
				int[] actuals = t.reverse(out);
				System.out.println(Arrays.toString(expecteds));
				System.out.println(Arrays.toString(actuals));
				assertArrayEquals(expecteds, actuals);
			}
		}
	}

}
