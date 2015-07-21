package org.rcsb.project3;

import java.io.Serializable;
import javax.vecmath.Point3d;
/**
 * This class generate fingerprint that is a sequence of angles which is the angle 
 * between each three points.
 * 
 * @author Chris Li
 */
public class StructuralAlphabetFingerprint implements SequenceFingerprint, Serializable {

	private static final long serialVersionUID = 1L;
	double[][] references = {
			{ 41.14,   75.53,  13.92,  -99.80,  131.88,  -96.27, 122.08,  -99.68},
			{108.24,  -90.12, 119.54,  -92.21,  -18.06, -128.93, 147.04,  -99.90},
			{-11.61, -105.66,  94.81, -106.09,  133.56, -106.93, 135.97, -100.63},
			{141.98, -112.79, 132.20, -114.79,  140.11, -111.05, 139.54, -103.16},
			{133.25, -112.37, 137.64, -108.13,  133.00,  -87.30, 120.54,   77.40},
			{116.40, -105.53, 129.32,  -96.68,  140.72,  -74.19, -26.65,  -94.51},
			{  0.40,  -81.83,   4.91, -100.59,   85.50,  -71.65, 130.78,   84.98},
			{119.14, -102.58, 130.83,  -67.91,  121.55,   76.25,  -2.95,  -90.88},
			{130.68,  -56.92, 119.26,   77.85,   10.42,  -99.43, 141.40,  -98.01},
			{114.32, -121.47, 118.14,   82.88, -150.05,  -83.81,  23.35,  -85.82},
			{117.16,  -95.41, 140.40,  -59.35,  -29.23,  -72.39, -25.08,  -76.16},
			{139.20,  -55.96, -32.70,  -68.51,  -26.09,  -74.44, -22.60,  -71.74},
			{-39.62,  -64.73, -39.52,  -65.54,  -38.88,  -66.89, -37.76,  -70.19},
			{-35.34,  -65.03, -38.12,  -66.34,  -29.51,  -89.10,  -2.91,   77.90},
			{-45.29,  -67.44, -27.72,  -87.27,    5.13,   77.49,  30.71,  -93.23},
			{-27.09,  -86.14,   0.30,   59.85,   21.51,  -96.30, 132.67,  -92.91}
		};
	String blockName = "abcdefghijklmnop";
	

	public StructuralAlphabetFingerprint() {
	}
	
	/**
     * Returns a fingerprint as sequence for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return
     */
	public StructuralAlphabetFeature getFingerprint(Point3d[] coordinates) {
		return null;
		
	}
	
	public String assign(double[] block) {
		double minRmsda = 9999999;
		char minBlock = 0;
		for (int i = 0; i < references.length; i++) {
			double rmsda = rmsda(references[i],block);
			if (rmsda < minRmsda) {
				minRmsda = rmsda;
				minBlock = blockName.charAt(i);
			}
		}
		return Character.toString(minBlock);
	}

	public double rmsda(double[] ref, double[] block) {
		double sum = 0;
		for (int i = 0; i < ref.length; i++) {
			double dif = angleModule360(ref[i] - block[i]);
			sum += dif * dif;
		}
		return sum;
	}
	
	public double angleModule360(double angle) {
		if (angle > 180)
			return angle - 360;
		else if (angle < -180)
			return angle + 360;
		else
			return angle;
	}
	
	public String getName() {
		return this.getClass().getSimpleName();
	}
}
