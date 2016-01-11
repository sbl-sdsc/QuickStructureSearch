package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.Function;

public class RadiusOfGyrationMapper implements Function<Point3d[], Float>{
	private static final long serialVersionUID = 1L;

	@Override
	public Float call(Point3d[] coordinates) throws Exception {
		double sum = 0;
		int n = 0;
		
	    for (int i = 0; i < coordinates.length-1; i++) {
	    	if (coordinates[i] == null) {
	    		continue;
	    	} else {
	    		n++;
	    	}
	    	for (int j = i; j < coordinates.length; j++) {
	    		if (coordinates[j] == null) continue;
	    		sum+= coordinates[i].distanceSquared(coordinates[j]);
	    	}
	    }

	    return (float)Math.sqrt(sum/n);
	}
}
