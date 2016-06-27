package org.rcsb.project3;

import java.util.List;
import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This is an interface for alignment algorithm that compare two proteins' fingerprint
 * 
 * @author Chris Li
 */
public interface AlignmentAlgorithmInterface extends PairFunction<Tuple2<String,String>,String,Float>{
	
	/**
	 * Returns the name of the algorithm
	 * @return algorihm name
	 */
	public String getName();
	
	/**
	 * Sets sequence(fingerprint) data
	 * @param sequences
	 */
	public void setSequence(Broadcast<Map<String,SequenceFeatureInterface<?>>> sequences);
	/**
	 * Sets coordinates
	 * @param coords
	 */
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> coords);
}
