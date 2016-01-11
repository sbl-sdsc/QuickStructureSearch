package org.rcsb.fingerprints;

import java.io.Serializable;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;


/**
 * This class generates a fingerprint (signature) for protein chains based upon
 * the writhingNumber Calculation method from Michael Levitt's paper
 *
 * @author Xueyang Li
 */
public class WrithingNumberFingerprint implements GenericFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	Vector3d rij1 = new Vector3d();
	Vector3d rij = new Vector3d();
	Vector3d ri1j1 = new Vector3d();
	Vector3d ri1j = new Vector3d();
	Vector3d a = new Vector3d();
	Vector3d b = new Vector3d();
	Vector3d c = new Vector3d();
	Vector3d d = new Vector3d();

    /**
     * Default constructor uses default parameters
     */
    public WrithingNumberFingerprint() {}
    
    public String getName() {
    	return this.getClass().getSimpleName();
    }
    
    public double[] getFingerprint(Point3d[] coords) {
    	// simplified version that avoids redundant calculation.
    	// can be further optimized by caching intermediate results
    	double sum = 0;
    	
    	for(int i = 0; i < coords.length - 2; i++)
    	{
    		if(coords[i] == null || coords[i+1] == null) continue;
    		// need to start at i + 2, otherwise coords[i+1] == coords[j] causes divide by zero
    		for(int j = i + 2; j < coords.length - 1; j++)
    		{
    			if(coords[j] == null || coords[j+1] == null) continue;
    			
    			// Vector a			
    			rij1.set(coords[i]);
    			rij1.sub(coords[j+1]);
    				
    			rij.set(coords[i]);
    			rij.sub(coords[j]);   			
    			
    			a.cross(rij1,rij);
    			a.normalize();
    			
    			// Vector b						
    			ri1j1.set(coords[i+1]);
    			ri1j1.sub(coords[j+1]);   			
    			
    			b.cross(rij1,ri1j1);
    			b.normalize();
    		
    			// Vector c i+1 can be j ??
    			ri1j.set(coords[i+1]);
    			ri1j.sub(coords[j]);
    			
    			c.cross(ri1j,ri1j1);
    			c.normalize();

    			// Vector d		    			
    			d.cross(ri1j,rij); 
    			d.normalize();
    			
    			// Get the dot product of each pair of vectors
    	
    			double ADot = a.dot(d);
    			double BDot = b.dot(a);
    			double CDot = c.dot(b);
    			double DDot = d.dot(c);
    			
    			// Calculate the value of angles:A,B,C,D
    			double A = Math.acos(ADot);
    			double B = Math.acos(BDot);
    			double C = Math.acos(CDot);
    			double D = Math.acos(DDot);
    			
    			// Calculate the Omega value
    			sum += A + B + C + D - 2*Math.PI;
    		}
    	}

    	double[] features = new double[1];
    	features[0] = sum/(2*Math.PI);

    	return features;
    }
    
}