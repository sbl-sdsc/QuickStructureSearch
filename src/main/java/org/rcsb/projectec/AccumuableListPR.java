package org.rcsb.projectec;

import javax.vecmath.Point3d;

import org.apache.spark.AccumulableParam;
import org.rcsb.structuralAlignment.SuperPosition;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import java.util.List;

public class AccumuableListPR implements AccumulableParam<List<Point3d[]>,Point3d[]>{
	private static final long serialVersionUID = 8776704701706611288L;
	private double rmsdThreshold;
	
	// create new SuperPositionQCP object that uses pre-centered coordinates
	private SuperPositionQCP qcp = new SuperPositionQCP(true); 
	
	public AccumuableListPR(double rmsdThreshold) {
		this.rmsdThreshold = rmsdThreshold;
	}

	@Override
	public List<Point3d[]> addAccumulator(List<Point3d[]> lib, Point3d[] fragment) {
		// pre-center fragment
		Point3d[] cFragment = SuperPosition.clonePoint3dArray(fragment);
		SuperPositionQCP.center(cFragment);
		
		// compare each archetype with the new fragment
		for (Point3d[] archetype: lib) {		
			qcp.set(archetype, cFragment);
			double rmsd = qcp.getRmsd();
			
			// if this fragment is similar to an existing fragment, don't add it to the library
			if (rmsd < rmsdThreshold) {
				return lib;
			}
		}
		
		// we found a new archetype, add it to the library
		lib.add(cFragment);
		System.out.println("Adding archetype: " + lib.size());
		return lib;
	}

	@Override
	public List<Point3d[]> addInPlace(List<Point3d[]> lib1, List<Point3d[]> lib2) {
		lib1.addAll(lib2);
		return lib1;
	}

	@Override
	public List<Point3d[]> zero(List<Point3d[]> lib) {
		lib.clear();
		return lib;
	}
}
