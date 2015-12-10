package org.rcsb.compress;

import java.io.Serializable;

public class CombinedTransform implements IntegerTransform, Serializable {
	private static final long serialVersionUID = 1L;
	private IntegerTransform transform1;
	private IntegerTransform transform2;
	
	public CombinedTransform(IntegerTransform transform1, IntegerTransform transform2) {
		this.transform1 = transform1;
		this.transform2 = transform2;
	}
	

	@Override
	public int[] forward(int[] data) {
		return transform2.forward(transform1.forward(data));
	}

	@Override
	public int[] reverse(int[] data) {
		return transform1.reverse(transform2.reverse(data));
	}
	
}
