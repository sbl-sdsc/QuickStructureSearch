package org.rcsb.projectm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.FlatMapFunction;
/**
 * Takes protein and returns a List of its fragments (which are Point3d arrays).
 * @author Emilia Copic
 *
 */
public class FlatMapToFragments implements FlatMapFunction<Point3d[], Point3d[]> {
	int fragmentSize = 6;

	public FlatMapToFragments(int fragSize) {
		this.fragmentSize = fragSize;
	}

	@Override
	public Iterable<Point3d[]> call(Point3d[] protein) throws Exception {
		List<Point3d[]> fragList = new ArrayList<>();
		for (int i = 0; i < protein.length - fragmentSize + 1; i++) {
			Point3d[] fragment = Arrays.copyOfRange(protein, i, i + fragmentSize);
			//get rid of nulls
			boolean result = true;
			for (int j = 0; j < fragment.length; j++) {
				result = result && (fragment[j] != null);
			}
			//if no nulls, add
			if (result) {
				fragList.add(fragment);
			}

		}
		return fragList;
	}

}
