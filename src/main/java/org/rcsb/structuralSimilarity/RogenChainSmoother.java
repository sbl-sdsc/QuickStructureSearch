package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;
import java.io.Serializable;

public class RogenChainSmoother implements ChainSmoother , Serializable
{
	private static final long serialVersionUID = 1L;
	private int iterations;
	public RogenChainSmoother(int iter){
		iterations = iter;
	}
	
	@Override
	public Point3d[] getSmoothedPoints(Point3d[] points) {
		// TODO Auto-generated method stub
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
		for (int i = start; i <= end; i++) {
			x[i] = points[i].x;
			y[i] = points[i].y;
			z[i] = points[i].z;
		}
		
		for (int i = 0; i < iterations; i++) {
			x = smoothRogen(x);
			y = smoothRogen(y);
			z = smoothRogen(z);
		}

		Point3d[] smoothedPoints = new Point3d[x.length];
		for (int i = 0; i < x.length; i++) {
			smoothedPoints[i] = new Point3d(x[i],y[i],z[i]);
		}
		
		return smoothedPoints;
	}

	/**
	 * P. Rogen, Evaluating protein structure descriptors and tuning Gauss integral
	 * based descriptors, J. Phys.: Condens. Matter (2005) 17, S1523-S1538
	 * @param x
	 * @return
	 */
	private static double[] smoothRogen(double[] x) {	    	
		double[] y = new double[x.length];
        double a = 2.4;
        double b = 2.1;
        
		y[0] = x[0];
		y[1] = (a*x[0] + b*x[1] + a*x[2])/(a+a+b);
		for (int i = 2; i < x.length-2; i++) {
			y[i] = (x[i-2]+ a*x[i-1] + b*x[i] + a*x[i+1] + x[i+2])/(2+a+a+b);
		}
		y[x.length-2] =  (a*x[x.length-3] + b*x[x.length-2] + a*x[x.length-1])/(a+a+b);
		y[x.length-1] = x[x.length-1];

		return y;
	}
}
