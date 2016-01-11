package org.rcsb.compress;

import java.io.Serializable;
import java.util.Arrays;

import me.lemire.integercompression.DeltaZigzagVariableByte;
import me.lemire.integercompression.IntWrapper;

public class PFORTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private DeltaZigzagVariableByte codec;
	

	@Override
	public int[] forward(int[] data) {
		int[] out = compress(data);
//		System.out.println("forward: " + data.length + " > " + out.length);
//		System.out.println(Arrays.toString(data));
//		System.out.println(Arrays.toString(out));
        return out;
	}

	@Override
	public int[] reverse(int[] data) {
		return uncompress(data);
	}
	
	private int[] compress(int[] data) {
		int[] outBuf = new int[data.length + 1024];
		IntWrapper inPos = new IntWrapper();
		IntWrapper outPos = new IntWrapper();
		codec.compress(data, inPos, data.length, outBuf, outPos);
		int[] out = new int[outPos.get()+1];
		System.arraycopy(outBuf, 0, out, 1, outPos.get());
		out[0] = data.length;
		return out;
	}
	
	private int[] uncompress(int[] data) {
		int [] outBuf = new int[data[0]+1024];
//		System.out.println("uncompress length: " + data.length);
		IntWrapper inPos = new IntWrapper(1);
		IntWrapper outPos = new IntWrapper(0);
		codec.uncompress(data, inPos, data.length-1, outBuf, outPos);
		return Arrays.copyOf(outBuf,  data[0]);
	}
	
}
