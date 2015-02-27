package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.rcsb.fingerprints.LinearFingerprint;

import scala.Tuple2;

public class ChainToLinearFeatureMapper implements PairFunction<Tuple2<String,Point3d[]>, String, Vector> {
	private static final long serialVersionUID = 1L;
	private LinearFingerprint fingerPrint;

	public ChainToLinearFeatureMapper(LinearFingerprint fingerPrint) {
		this.fingerPrint = fingerPrint;
	}

	@Override
	public Tuple2<String, Vector> call(Tuple2<String,Point3d[]> tuple) throws Exception {
		int[] intFeatures = this.fingerPrint.getLinearFingerprint(tuple._2);
		double[] features = new double[intFeatures.length];
		for (int i = 0; i < features.length; i++) {
			features[i] = intFeatures[i];
		}
		return new Tuple2<String, Vector>(tuple._1, Vectors.dense(features));
	}
}
