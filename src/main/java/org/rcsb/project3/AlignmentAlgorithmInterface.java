package org.rcsb.project3;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This is an interface for alignment algorithm that compare two proteins' fingerprint
 * 
 * @author Chris Li
 */
public interface AlignmentAlgorithmInterface extends PairFunction<Tuple2<Integer,Integer>,String,Float>{
	/**
	 * load seqeuence(fingerprint) data
	 * @param sequences
	 */
	public void setSequence(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> sequences);
	/**
	 * load coordinates
	 * @param coords
	 */
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> coords);
}
