package org.rcsb.project4;

import org.apache.spark.AccumulatorParam;

/**
 * This class implements the AccumulatorParam of Spark for Object type of Long
 * It is used to record the time used in parallel threads
 * 
 * @author Chris Li
 */
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
