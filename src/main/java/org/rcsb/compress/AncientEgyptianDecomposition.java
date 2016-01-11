package org.rcsb.compress;

import java.io.Serializable;

/**
 * Code obtained from https://github.com/cscheiblich/JWave and adopted to
 * new infrastructure.
 * 
 * JWave is distributed under the MIT License (MIT); this file is part of.
 *
 * Copyright (c) 2008-2015 Christian Scheiblich (cscheiblich@gmail.com)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * A wavelet transform method for arrays and signals of arbitrary lengths, even
 * odd lengths. The array is decomposed in several parts of optimal lengths by
 * applying the ancient Egyptian decomposition. Hereby, the array or signal is
 * decomposed to the largest possible sub arrays of two the power of p.
 * Afterwards each sub array is transformed forward and copied back to the
 * discrete position of the input array. The reverse transform applies the same
 * vice versa. In more detail the ancient Egyptian Multiplication can be easily
 * explained by the following example: 42 = 2^5 + 2^3 + 2^1 = 32 + 8 + 2.
 * However, an array or signal of odd length produces the smallest ancient
 * Egyptian multiplier 2^0 which is actually 1. Therefore, the matching sub
 * array or signal is untouched an the coefficient is actually the wavelet
 * coefficient of wavelet space of level 0. For an "orthonormal" wavelet this
 * holds. See: http://en.wikipedia.org/wiki/Ancient_Egyptian_multiplication
 * 
 * @date 14.08.2010 10:43:28
 * @author Christian Scheiblich (cscheiblich@gmail.com)
 * @author Peter Rose (adopted to IntegerTransform)
 */
public class AncientEgyptianDecomposition implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private IntegerTransform transform;

	public AncientEgyptianDecomposition(IntegerTransform transform) {
		this.transform = transform;
	}

	/**
	 * This forward method decomposes the given array of arbitrary length to sub
	 * arrays while applying the ancient Egyptian decomposition. Each sub array is
	 * transformed by the selected basic transform and the resulting wavelet
	 * coefficients are copied back to their original discrete positions.
	 * 
	 * @date 14.08.2010 10:43:28
	 * @author Christian Scheiblich (cscheiblich@gmail.com)
	 * @author Peter Rose (adopted for IntegerTransformation)
	 */
	public int[] forward(int[] in) {

		int[] ancientEgyptianMultipliers = MathToolKit.decompose(in.length);

		int[] out = new int[in.length];
		int offSet = 0;

		for (int m: ancientEgyptianMultipliers) {
			
			int len = (int)MathToolKit.scalb(1., m);
			int[] tmp = new int[len];

			for(int i = 0; i < len; i++) {
				tmp[i] = in[i + offSet];
			}

			tmp = transform.forward(tmp);

			for (int i = 0; i < len; i++) {
				out[i + offSet] = tmp[i];
			}

			offSet += tmp.length;

		}

		return out;

	}

	/**
	 * This reverse method awaits an array of arbitrary length in wavelet space
	 * keeping the wavelet already decomposed by the ancient Egyptian
	 * decomposition. Therefore, each of the existing sub arrays of length 2^p is
	 * reverse transformed by the selected basic transform and the resulting
	 * coefficients of time domain are copied back to their original discrete
	 * positions.
	 * 
	 * @date 14.08.2010 10:43:28
	 * @author Christian Scheiblich (cscheiblich@gmail.com)
	 * @author Peter Rose (adopted to IntegerTransform)
	 */
	public int[] reverse(int[] in) {
		
		int[] ancientEgyptianMultipliers = MathToolKit.decompose(in.length);
		
		int[] out = new int[in.length];
		int offSet = 0;

		for (int m: ancientEgyptianMultipliers) {
			int len = (int)MathToolKit.scalb(1., m);

			int [] tmp = new int[len];

			for (int i = 0; i < len; i++ ) {
				tmp[i] = in[i + offSet];
			}

			tmp = transform.reverse(tmp);

			for (int i = 0; i < len; i++) {
				out[i + offSet] = tmp[i];
			}

			offSet += tmp.length;

		}

		return out;

	}
}
