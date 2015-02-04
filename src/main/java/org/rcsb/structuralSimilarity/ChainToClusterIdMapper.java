package org.rcsb.structuralSimilarity;

import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * Maps a chains represented by chainId and coordinates to a sequence cluster id.
 * It use the cluster map passed into the constructor to determine the cluster id.
 * 
 * @author Peter Rose
 *
 */
public class ChainToClusterIdMapper implements PairFunction<Tuple2<String,Integer>,Integer, Tuple2<String,Integer>> {
	private static final long serialVersionUID = 1L;
	private Broadcast<Map<String,Integer>> clusterMap;

	public ChainToClusterIdMapper(Broadcast<Map<String,Integer>> clusterMap) {
		this.clusterMap = clusterMap;
	}

	/**
	 * Maps a <String, Integer> to a <clusterIndex, <String, Integer>> pair
	 */
	@Override
	public Tuple2<Integer, Tuple2<String, Integer>> call(
			Tuple2<String, Integer> t) throws Exception {
		Integer clusterId = clusterMap.getValue().get(t._1);
		if (clusterId == null) {
			clusterId = -1;
		}
		return new Tuple2<Integer, Tuple2<String, Integer>>(clusterId, t);
	}
}
