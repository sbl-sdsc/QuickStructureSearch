package org.rcsb.structuralAlignment;

import javax.vecmath.Point3d;

/**
 * This class calculates the distance RMSD between two Point3d arrays.
 * 
 * @author Dane Malangone
 * @author Justin Li
 * @author Reyd Nguyen
 * @author Joe Sun
 */
public class DistanceRmsd {
	
	/**
	 * Returns the distance RMSD between two Point3d arrays
	 * 
	 * @param x the first Point3d array
	 * @param y the second Point3d array
	 * @return distance RMSD
	 */
	public static double getDistanceRmsd (Point3d[] x, Point3d[] y) {

		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y length needs to be identical. x.length = " + x.length + " and y.length = " + y.length);	
		}
		double drmsd = 0.0;
		double sum = 0.0;
		int length = x.length;
		for (int i = 0; i < length-1; i++) {
			for (int j = i + 1; j < length; j++) {
				double dist1 = x[i].distance(x[j]);
				double dist2 = y[i].distance(y[j]);
				sum += Math.pow((dist1 - dist2),2);
			}
		}
		drmsd = sum * 2 / (length * (length-1));
		drmsd = Math.sqrt(drmsd);

		return drmsd;
	}
}
