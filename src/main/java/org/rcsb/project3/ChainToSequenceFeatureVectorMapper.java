package org.rcsb.project3;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

/**
 * This class is used to map the chain to sequence feature
 * @author Chris Li
 *
 */
public class ChainToSequenceFeatureVectorMapper implements PairFunction<Tuple2<String,Point3d[]>, String, SequenceFeatureInterface<?>> {
	private static final long serialVersionUID = 1L;
	private SequenceFingerprint fingerPrint;

	public ChainToSequenceFeatureVectorMapper(SequenceFingerprint fingerPrint) {
		this.fingerPrint = fingerPrint;
	}

	@Override
	public Tuple2<String, SequenceFeatureInterface<?>> call(Tuple2<String,Point3d[]> tuple) throws Exception {
		SequenceFeatureInterface<?> features = this.fingerPrint.getFingerprint(tuple._2);
		return new Tuple2<String, SequenceFeatureInterface<?>>(tuple._1, features);
	}
}
