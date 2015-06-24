package org.rcsb.project4;

import org.apache.spark.AccumulatorParam;

public class TimeAccumulator implements AccumulatorParam<Long>{

	private static final long serialVersionUID = 1L;

	@Override
	public Long addInPlace(Long r, Long t) {
		return r + t;
	}

	@Override
	public Long zero(Long arg0) {
		return 0L;
	}

	@Override
	public Long addAccumulator(Long r, Long t) {
		return r + t;
	}

}
