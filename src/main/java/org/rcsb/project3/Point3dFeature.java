package org.rcsb.project3;

import java.io.Serializable;

import javax.vecmath.Point3d;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for comparing protein's coordinates directly.
 *  
 * @author Chris Li
 */
public class Point3dFeature implements SequenceFeatureInterface<Point3d>, Serializable {

	private static final long serialVersionUID = 1L;
	private Point3d[] coords;
	// Some setting for calculate similarity score 
	private double match = 1;
	private double mismatch = -1;
	private double gap = -1;
	private double diff = 4;
	
	public Point3dFeature(Point3d[] coords) {
		this.coords = coords;
	}
	
	@Override
	public double similarity(SequenceFeatureInterface<Point3d> sequence2,
			int i, int j) {
		// check for gap
		if (this.get(i) == null || sequence2.get(j) == null)
			return gap;
		else {
			double dist = this.get(i).distance(sequence2.get(j));
			if (dist < diff)
				return match;
			else {
				double score = 1 - dist/(diff*2.5);
				if (score > 0)
					return score;
				else
					return mismatch;
			}
		}
	}

	@Override
	public boolean identity(SequenceFeatureInterface<Point3d> sequence2, int i,
			int j) {
		if (this.get(i).distance(sequence2.get(j)) < diff)
			return true;
		else
			return false;
	}

	@Override
	public Point3d[] getSequence() {
		return coords;
	}

	@Override
	public Point3d get(int index) {
		return coords[index];
	}

	@Override
	public int length() {
		return coords.length;
	}

	@Override
	public String toString(int index) {
		return null;
	}

	@Override
	public double todouble(int index) {
		return 0;
	}

	@Override
	public Point3d[] getCoords() {
		return coords;
	}

}
