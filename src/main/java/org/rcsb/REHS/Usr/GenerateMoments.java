package org.rcsb.REHS.Usr;

import javax.vecmath.Point3d;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * Generate the USR moments for a given {@link Point3d} array.
 * @author Anthony Bradley, Michael Wang
 *
 */
public class GenerateMoments {
	


	/**
	 * Generate the USR moments for a given input array of 3D points.
	 * @param inputArray the array of {@link Point3d}
	 * @return the moments of this array
	 */
//	private static int numPoints = 6;
	public static double[] getMoments(Point3d[] inputArray) {
		double[] outArray = new double[12];
		Point3d[] Points = getFourPoints(inputArray);
		for (int i=0; i<4; i++) {
			float[] threeMoments = getThreeMoments(getDistribution(Points[i], inputArray));
			for (int j=0; j<3; j++){
				outArray[i*3+j] = threeMoments[j];
			}
		}
		return outArray;
	}

	/**
	 * Get the four points of interest. The middle point, the point furthes
	 * @param inputArray the inputArray of points
	 * @return
	 */

	private static Point3d[] getFourPoints(Point3d[] inputArray) {		
		Point3d[] outArray = new Point3d[4];
		outArray[0] = getCentroid(inputArray);
		outArray[1] = getFarthestPoint(inputArray, outArray[0], null);
		for (int i=2; i<4; i++){
			outArray[i] = getFarthestPoint(inputArray, outArray[i-1], outArray[i-2]);
		}
		return outArray;
	}

	/**
	 * Function to get the farthest point from a single point in an array
	 * of points.
	 * @param inputArray the input {@link Point3d} array
	 * @param queryPoint the point to be farthest from
	 * @return the farthest point
	 */
	private static Point3d getFarthestPoint(Point3d[] inputArray, Point3d queryPoint, Point3d reUsed) {
		double maxDist = -1.0;
		Point3d maxPoint = null;
		for(Point3d point3d : inputArray) {
			if(point3d != null){
				double currentDist = point3d.distance(queryPoint);
				if(point3d != reUsed && currentDist>maxDist){
					maxPoint = point3d;
					maxDist = currentDist;  
			}
			}
		}
		return maxPoint;
	}

	/**
	 * Get the centroid of a list of {@link Point3d} objects.
	 * @param points the input {@link Point3d} objects
	 * @return the centroid
	 */
	private static Point3d getCentroid(Point3d[] points) {
		Point3d centroid = new Point3d();
		double sumX = 0;
		double sumY = 0;
		double sumZ = 0;
		final int nPoints = points.length;
		for (int n = 0; n < nPoints; n++) {
			
			if(points[n] != null)
			{
			sumX += (double) points[n].x;
			sumY += (double) points[n].y;
			sumZ += (double) points[n].z;
			}
		}
		centroid.x = sumX / nPoints;
		centroid.y = sumY / nPoints;
		centroid.z = sumZ / nPoints;

		return centroid;
	}
	/**
	 * Get the mean, variance and the third central moment from an array of doubles.
	 * @param distribution the array of doubles
	 * @return a float of length three. Encodes the mean, the variance and the
	 * third central moment
	 */
	private static float[] getThreeMoments(double[] distribution) {
		DescriptiveStatistics descriptiveStatistics = new DescriptiveStatistics();
		for (double val : distribution){
			if(val != 0)
				descriptiveStatistics.addValue(val);
		}
		float[] outVals = new float[3];
		outVals[0] = (float) descriptiveStatistics.getMean();
		outVals[1] = (float) descriptiveStatistics.getVariance();
		outVals[2] = (float) descriptiveStatistics.getSkewness();
		return outVals;
	}

	/**
	 * Get the distribution of euclidean distances between a point 
	 * and a series of points.
	 * NB Currently calculates all (including itself) so there will always be a 
	 * zero.
	 * @param singlePoint3d the input point
	 * @param inputArray the points to find distances to
	 * @return the distribution of distances as an array of floats
	 */
	private static double[] getDistribution(Point3d singlePoint3d, Point3d[] inputArray) {
		double[] outArray = new double[inputArray.length];
		for(int i=0; i<inputArray.length;i++){
			if(inputArray[i] != null)
			{
			outArray[i] = inputArray[i].distance(singlePoint3d);
			}
		}
		return outArray;
	}
}