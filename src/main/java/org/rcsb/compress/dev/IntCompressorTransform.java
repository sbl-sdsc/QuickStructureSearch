package org.rcsb.compress.dev;

import java.io.Serializable;

import org.rcsb.compress.IntegerTransform;

import me.lemire.integercompression.IntCompressor;
import me.lemire.integercompression.SkippableIntegerCODEC;

public class IntCompressorTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private IntCompressor codec;
	
	public IntCompressorTransform(SkippableIntegerCODEC codec) {
		this.codec = new IntCompressor(codec);
	}
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public int[] forward(int[] data) {
		return codec.compress(data);
	}

	@Override
	public int[] reverse(int[] data) {
		return codec.uncompress(data);
	}
	
}
