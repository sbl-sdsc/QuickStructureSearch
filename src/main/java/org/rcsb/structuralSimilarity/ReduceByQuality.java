package org.rcsb.structuralSimilarity;

import org.apache.spark.api.java.function.Function2;

import scala.Tuple2;

public class ReduceByQuality implements Function2<Tuple2<String, Integer>, Tuple2<String, Integer>, Tuple2<String, Integer>> {
	private static final long serialVersionUID = 1L;

	@Override
	public Tuple2<String, Integer> call(Tuple2<String, Integer> v1,
			Tuple2<String, Integer> v2) throws Exception {
		return v1._2 >= v2._2? v1 : v2;
	}

}
