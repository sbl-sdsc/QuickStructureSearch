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
			}
		};

		return Math.sqrt(sum * 2.0 / (x.length * (x.length-1)));
	}
	
	
	/**
	 * Returns and approximate distance RMSD between two Point3d arrays
	 * Fit with adRMSD = 1.3325 * x - 0.0856, R^2 = 0.961
	 * @param x the first Point3d array
	 * @param y the second Point3d array
	 * @return distance RMSD
	 */
	public static double getDistanceRmsdApproximation(Point3d[] x, Point3d[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y length needs to be identical. x.length = " + x.length + " and y.length = " + y.length);	
		}

		double sum = 0.0;
		
	    double diff = x[0].distance(x[7]) - y[0].distance(y[7]);
		sum += diff * diff;
		diff = x[0].distance(x[3]) - y[0].distance(y[3]);
		sum += diff * diff;
		diff = x[1].distance(x[6]) - y[1].distance(y[6]);
		sum += diff * diff;
		diff = x[2].distance(x[5]) - y[2].distance(y[5]);
		sum += diff * diff;
		diff = x[4].distance(x[7]) - y[4].distance(y[7]);
		sum += diff * diff;
		diff = x[1].distance(x[5]) - y[1].distance(y[5]);
		sum += diff * diff;
		diff = x[2].distance(x[6]) - y[2].distance(y[6]);
		sum += diff * diff;
		
//		return (Math.sqrt(sum/7) + 0.0856)/1.3325;
//		return Math.sqrt(sum/7) * 0.75;
		return Math.sqrt(sum * 0.08357);
	}
	
	public static double getDistanceRmsdScaled(Point3d[] x, Point3d[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y length needs to be identical. x.length = " + x.length + " and y.length = " + y.length);	
		}

		double sum = 0.0;
		
		for (int i = 0; i < x.length-1; i++) {
			for (int j = i + 1; j < x.length; j++) {
				double diff = x[i].distance(x[j]) - y[i].distance(y[j]);
				sum += diff * diff;
			}
		}

		double drmsd = sum * 2 / (x.length * (x.length-1));
		drmsd = Math.sqrt(drmsd);
		
		double xl = Math.sqrt(x[0].distance(x[x.length-1]));
		double yl = Math.sqrt(y[0].distance(y[x.length-1]));
		double ml = (xl + yl) * 0.5;
		drmsd *= ml;

		return drmsd;
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
