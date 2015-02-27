package org.rcsb.structuralAlignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.PairFlatMapFunction;

import scala.Tuple2;

public class ChainToFragmentMapper implements PairFlatMapFunction<Tuple2<String, Point3d[]>, String, Point3d[]>{
	private static final long serialVersionUID = 1L;
	private int fragmentLength;

	public ChainToFragmentMapper(int fragmentLength) {
		this.fragmentLength = fragmentLength;
	}
	

	@Override
	public Iterable<Tuple2<String, Point3d[]>> call(Tuple2<String, Point3d[]> tuple) {
		String chainId = tuple._1;
		Point3d[] points = tuple._2;
		List<Tuple2<String, Point3d[]>>list = new ArrayList<>(points.length-this.fragmentLength);
		
		for (int i = 0; i < points.length-this.fragmentLength; i++){
    		if (hasGaps(points, i)){
    			continue; //Skip if gap exists
    		}
    		Point3d[] fragment = Arrays.copyOfRange(points, i, i+this.fragmentLength);
    		list.add(new Tuple2<String, Point3d[]>(chainId, fragment));
		}
		return list;
	}

	/**
	 * Returns true if there is a gap between the C alpha atoms
	 * within a fragment. Note, the first position is not checked, 
	 * since the gap info is for a gap between residue i and i + 1.
	 * @param gaps true if there is a gap in the fragment
	 * @param index start residue of fragment
	 * @return
	 */
	private boolean hasGaps(Point3d[] coords, int index) {
		for (int i = index; i < index+this.fragmentLength; i++) {
			if (coords[i] == null) {
				return true;
			}
		}
		return false;
	}
}
