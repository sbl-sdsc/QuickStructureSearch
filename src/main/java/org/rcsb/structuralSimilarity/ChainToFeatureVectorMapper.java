package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.rcsb.fingerprints.GenericFingerprint;

import scala.Tuple2;

public class ChainToFeatureVectorMapper implements PairFunction<Tuple2<String,Point3d[]>, String, Vector> {
	private static final long serialVersionUID = 1L;
	private GenericFingerprint fingerPrint;

	public ChainToFeatureVectorMapper(GenericFingerprint fingerPrint) {
		this.fingerPrint = fingerPrint;
	}

	@Override
	public Tuple2<String, Vector> call(Tuple2<String,Point3d[]> tuple) throws Exception {
		double[] features = this.fingerPrint.getFingerprint(tuple._2);		
		return new Tuple2<String, Vector>(tuple._1, Vectors.dense(features));
	}
}
