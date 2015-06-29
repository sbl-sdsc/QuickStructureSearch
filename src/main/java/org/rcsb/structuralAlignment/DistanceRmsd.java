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

		double sum = 0.0;
		
		for (int i = 0; i < x.length-1; i++) {
			for (int j = i + 1; j < x.length; j++) {
				// faster version
				double diff = x[i].distance(x[j]) - y[i].distance(y[j]);
				sum += diff * diff;
				// original version
				//				double dist1 = x[i].distance(x[j]);
				//				double dist2 = y[i].distance(y[j]);
				//				sum += Math.pow((dist1 - dist2),2);
			}
		}

//		drmsd = sum * 2 / (length * (length-1));
//		drmsd = Math.sqrt(drmsd);

		return Math.sqrt(sum * 2.0 / (x.length * (x.length-1)));
	}

	/**
	 * Returns the centroid distance RMSD between two Point3d arrays
	 * @param x the first Point3d array
	 * @param y the second Point3d array
	 * @return centroid distance RMSD
	 */
	public static double getCentroidDistanceRmsd(Point3d[] x, Point3d[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y length needs to be identical. x.length = " + x.length + " and y.length = " + y.length);	
		}

		Point3d xCentroid = centroid(x);
		Point3d yCentroid = centroid(y);

		double sum = 0;
		for (int i = 0; i < x.length; i++) {
			double d = xCentroid.distance(x[i]) - yCentroid.distance(y[i]);
			sum += d * d;
		}

		return Math.sqrt(sum / x.length);
	}

	public static Point3d centroid(Point3d[] x) {
		Point3d center = new Point3d();
		for (Point3d p: x) {
			center.add(p);
		}
		center.scale(1.0/x.length);
		return center;
	}
}
