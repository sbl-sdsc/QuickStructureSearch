package org.rcsb.fingerprints;

import javax.vecmath.Point3d;

import java.io.Serializable;

public class PointToPointDistanceFingerprint implements GenericFingerprint, Serializable {
	
	private static final long serialVersionUID = 1L;
	private double binSize = 2.0;
	private double maxDist = 1000;
	private double minDist = 3.3;
	private int featureCount = (int)((this.maxDist-this.minDist)/this.binSize);
	
	/**
     * Default constructor uses default parameters
     */
    public PointToPointDistanceFingerprint() {}
    
	/**
     * Constructor with all parameters
     * @param max min fragment length
     */
    public PointToPointDistanceFingerprint (double max, double min, double binSize) {
        this.maxDist = max;
        this.minDist = min;
        this.binSize = binSize;
        this.featureCount = (int)((this.maxDist-this.minDist)/this.binSize);
	}
    
	@Override
	public double[] getFingerprint(Point3d[] coordinates) {
		// TODO Auto-generated method stub
		double[] features = new double[this.featureCount];
    	double scale = 1/this.binSize;
    	if (coordinates.length <= 1) {
    		return features;
    	}
    	for (int i = 0; i < coordinates.length-1; i++) {
    		Point3d p1 = coordinates[i];
    		// skip gaps
    		if (p1 == null)
    			continue;
    		for (int j = i+1; j <= coordinates.length-1; j++) {
    			Point3d p2 = coordinates[j];
    			// skip gaps
    			if (p2  == null) {
    				continue;
    			}
	    		// calculate point to point distance
	    		// and bin values
    			double dist = p1.distance(p2);
    			if (dist < minDist || dist > maxDist)
    				continue;
	    		int bin = (int)Math.round(scale*dist);
	    		if (bin > this.featureCount-1) {
	    			continue;
	    		}
	    		features[bin]++;
    		}
    	}
		return features;
	}

	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return this.getClass().getSimpleName();
	}


}
