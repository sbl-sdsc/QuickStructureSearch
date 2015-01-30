package org.rcsb.structuralSimilarity;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class implements a filter in a Spark workflow to filter out protein chains
 * that are not contained in the passed in chainId list.
 * 
 * @author Peter Rose
 */
public class ChainIdFilter<E> implements Function<Tuple2<String,E>, Boolean> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<String>> chainIds;

	/**
	 * @param list of chain ids
	 */
	public ChainIdFilter(Broadcast<List<String>> chainIds) {
		this.chainIds = chainIds;
	}

	@Override
	public Boolean call(Tuple2<String,E> tuple) throws Exception {
		Collection<String> values = chainIds.getValue();
		String s = tuple._1;
		
		// pdbId.chainId must contain a period
		if (!s.contains(".")) {
			return false;
		}

		if (s.contains(",")) {
			System.out.println(s);
			// comma separated list of chainIds
			String[] ids = s.split(",");
			for (String id: ids) {
				System.out.println(Arrays.toString(ids));
				if (!values.contains(id)) {
					return false;
				}
			}
			return true;
		} else {
			// single chainId
			return values.contains(tuple._1);
		}	
	}
}
