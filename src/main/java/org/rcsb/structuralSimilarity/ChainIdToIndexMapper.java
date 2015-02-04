package org.rcsb.structuralSimilarity;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
* This class maps a tuple of chainIds <String, String> to a tuple of chain indices <Integer, Integer>
* based on the passed in chainId list.
* 
* @author Peter Rose
*/
public class ChainIdToIndexMapper implements PairFunction<Tuple2<String, String>,Integer, Integer> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<String>> chainIds;
	
	public ChainIdToIndexMapper(Broadcast<List<String>> chainIds) {
		this.chainIds = chainIds;
	}

	/**
	 * Maps chainId pairs <String, String> to chainId indices
	 */
	@Override
	public Tuple2<Integer, Integer> call(Tuple2<String, String> tuple) throws Exception {	
		long index1 = chainIds.getValue().indexOf(tuple._1);
		long index2 = chainIds.getValue().indexOf(tuple._2);
		return new Tuple2<Integer,Integer>((int)index1, (int)index2);
	}
}
