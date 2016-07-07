package org.rcsb.projectva;
import java.util.List;

import javax.vecmath.Point3d;

import org.rcsb.project3.*;

import org.apache.spark.AccumulableParam;

public class ListAccumulator implements AccumulableParam<List<Point3d[]>, Point3d[]> {

	@Override
	public List<Point3d[]> addAccumulator(List<Point3d[]> arg0, Point3d[] arg1) {
		// TODO Auto-generated method stub
		arg0.add(arg1);
		return arg0;
	}

	@Override
	public List<Point3d[]> addInPlace(List<Point3d[]> arg0, List<Point3d[]> arg1) {
		// TODO Auto-generated method stub
		
		for (Point3d[] iter: arg1) {
			arg0.add(iter);
		}
		return arg0;
	}

	@Override
	public List<Point3d[]> zero(List<Point3d[]> arg0) {
		// TODO Auto-generated method stub
		
		return arg0;
	}

	

}
