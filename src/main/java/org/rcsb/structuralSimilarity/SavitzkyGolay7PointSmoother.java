package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import java.io.Serializable;

public class SavitzkyGolay7PointSmoother implements ChainSmoother , Serializable{
	private static final long serialVersionUID = 1L;
	private int iterations;
	private int[] coefficients = {-2,3,6,7,6,3,-2};
	private int norm = 21;
	public SavitzkyGolay7PointSmoother(int iter){
		iterations = iter;
	}

	@Override
	public Point3d[] getSmoothedPoints(Point3d[] points) {	
		double[] x = new double[points.length];
		double[] y = new double[points.length];
		double[] z = new double[points.length];

		int start = 0;
		int end = points.length-1;

		
		// skip N-terminal gap (start of chain)
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				start = i;
				break;
			}
		}
		
		// skip C-terminal gap (end of chain)
		for (int i = points.length-1; i > start; i--) {
			if (points[i] != null) {
				end = i;
				break;
			}
		}		
		
		// Only take in the points without front and end
		for (int i = start; i < end; i++) {
			x[i] = points[i].x;
			y[i] = points[i].y;
			z[i] = points[i].z;
		}
		
		for (int i = 0; i < iterations; i++) {
			x = smoothSavitzkyGolay(x, coefficients, norm);
			y = smoothSavitzkyGolay(y, coefficients, norm);
			z = smoothSavitzkyGolay(z, coefficients, norm);
		}

		// there will be fewer points after smoothing
		Point3d[] smoothedPoints = new Point3d[x.length];
		for (int i = 0; i < x.length; i++) {
			smoothedPoints[i] = new Point3d(x[i],y[i],z[i]);
		}
		
		return smoothedPoints;
	}
	
    private static double[] smoothSavitzkyGolay(double[] x, int[] coefficients, int norm) {
    	int n = coefficients.length - 1;
    	int m = x.length - n;
    	
    	double[] y = new double[m];
    	double[] p = new double[coefficients.length];
    	
        for (int i = 1; i < coefficients.length; i++) {
        	p[i] = x[i-1];	
        }
        for (int i = 0; i < m; i++) {
        	for (int k = 0; k < n; k++) {
        		p[k]= p[k+1];
        	}
        	p[n] = x[i+n];
        	double sum = 0;
        	for (int k = 0; k < coefficients.length; k++) {
        		sum += coefficients[k] * p[k];
        	}
        	y[i] = sum/norm;
        }
        return y;
    }
}
