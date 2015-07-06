package org.rcsb.project5;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.hadoop.io.SimplePolymerChain;

import scala.Tuple2;

/**
 * This class returns structural clusters given a sequence cluster.
 * 
 * @author Justin Li
 * @author Joe Sun
 * @author Peter Rose
 *
 */
public class StructureClusterer implements VoidFunction<Tuple2<Integer, Iterable<String>>> {
	private static final long serialVersionUID = 1L;
	private Broadcast<Map<String, SimplePolymerChain>> chainMap;
	private double maxRmsd;

	public StructureClusterer(Broadcast<Map<String, SimplePolymerChain>> chainMap, double maxRmsd) {
		this.chainMap = chainMap;
		this.maxRmsd = maxRmsd;
	}

	/**
	 * This method is currently not used
	 */
	@Override
	public void call(Tuple2<Integer, Iterable<String>> tuple) throws Exception {
		Map<String, SimplePolymerChain> map = chainMap.getValue();
//		List<Tuple2<String, Integer[]>> writerList = new ArrayList<Tuple2<String, Integer[]>>();
		int count = 1;

		System.out.println("*** Cluster: " + tuple._1 + ": ");
		StructuralClusterCreator clusterCreator = new StructuralClusterCreator();
		List<Tuple2<String, SimplePolymerChain>> list = new ArrayList<Tuple2<String, SimplePolymerChain>>();
		for (String id: tuple._2) {
			SimplePolymerChain c = map.get(id);
			if(c != null)
			{
				list.add(new Tuple2<String, SimplePolymerChain>(id, c));
			}
		}
		for (List<Tuple2<String, SimplePolymerChain>> a: clusterCreator.createStructuralCluster(list, 1.0)) {
			System.out.println("***** " + tuple._1 + "-" + count + " Structural Cluster: ");
			for(Tuple2<String, SimplePolymerChain> tuple2: a) {
				System.out.println(tuple2._1 + ": " + tuple2._2.getSequence());
//				writerList.add(new Tuple2<String, Integer[]>(tuple2._1, new Integer[]{tuple._1, count}));
			}
//			writeToCsv(writer, writerList);
			System.out.println();
			count++;
		}
		/*for (String id: tuple._2) {
			SimplePolymerChain c = map.get(id);
			System.out.println(id + ": " + c.getSequence());
		}*/
	}
	
	/**
	 * Returns a list of structural clusters given a sequence cluster
	 * @param tuple
	 * @return writerList
	 * @throws Exception
	 */
	public List<Tuple2<String, Integer[]>> getStructuralClusters(Tuple2<Integer, Iterable<String>> tuple) throws Exception {
		Map<String, SimplePolymerChain> map = chainMap.getValue();
		List<Tuple2<String, Integer[]>> writerList = new ArrayList<Tuple2<String, Integer[]>>();
		int count = 1;

		System.out.println("*** Cluster: " + tuple._1 + ": ");
		StructuralClusterCreator clusterCreator = new StructuralClusterCreator();
		List<Tuple2<String, SimplePolymerChain>> list = new ArrayList<Tuple2<String, SimplePolymerChain>>();
		for (String id: tuple._2) {
			SimplePolymerChain c = map.get(id);
			if(c != null)
			{
				list.add(new Tuple2<String, SimplePolymerChain>(id, c));
			} else{
				//adding the proteins with gaps
				writerList.add(new Tuple2<String, Integer[]>(id, new Integer[]{tuple._1, null}));
			}
		}
		for (List<Tuple2<String, SimplePolymerChain>> a: clusterCreator.createStructuralCluster(list, maxRmsd)) {
//			System.out.println("***** " + tuple._1 + "-" + count + " Structural Cluster: ");
			for(Tuple2<String, SimplePolymerChain> tuple2: a) {
//				System.out.println(tuple2._1 + ": " + tuple2._2.getSequence());
				writerList.add(new Tuple2<String, Integer[]>(tuple2._1, new Integer[]{tuple._1, count}));
			}
//			System.out.println();
			count++;
		}
		return writerList;
	}
}
