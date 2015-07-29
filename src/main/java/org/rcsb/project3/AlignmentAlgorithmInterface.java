package org.rcsb.project3;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

public interface AlignmentAlgorithmInterface extends PairFunction<Tuple2<Integer,Integer>,String,Float>{
	/**
	 * load seqeuence data
	 * @param data
	 */
	public void setSequence(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data);

	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> sequence);
}
