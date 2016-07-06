package org.rcsb.projectec;

import javax.vecmath.Point3d;

import org.apache.spark.AccumulableParam;

import java.util.List;


public class AccumuableList implements AccumulableParam<List<Point3d[]>,Point3d[]>{

	@Override
	public List<Point3d[]> addAccumulator(List<Point3d[]> lib, Point3d[] fragment) {
		lib.add(fragment);
		return lib;
	}

	@Override
	public List<Point3d[]> addInPlace(List<Point3d[]> arg0, List<Point3d[]> arg1) {
		// TODO Auto-generated method stub
		arg0.addAll(arg1);
		return arg0;
	}

	@Override
	public List<Point3d[]> zero(List<Point3d[]> arg0) {
		// TODO Auto-generated method stub
		return arg0;
	}

}
