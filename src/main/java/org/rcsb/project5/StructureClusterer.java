package org.rcsb.project5;

import java.util.Map;

import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.hadoop.io.SimplePolymerChain;

import scala.Tuple2;

public class StructureClusterer implements VoidFunction<Tuple2<Integer, Iterable<String>>> {
	private static final long serialVersionUID = 1L;
	private Broadcast<Map<String, SimplePolymerChain>> chainMap;

	public StructureClusterer(Broadcast<Map<String, SimplePolymerChain>> chainMap) {
		this.chainMap = chainMap;	
	}
	
	@Override
	public void call(Tuple2<Integer, Iterable<String>> tuple) throws Exception {
		Map<String, SimplePolymerChain> map = chainMap.getValue();
		
		System.out.println("*** Cluster: " + tuple._1 + ": ");

		for (String id: tuple._2) {
			SimplePolymerChain c = map.get(id);
			System.out.println(id + ": " + c);
		}
	}
}
