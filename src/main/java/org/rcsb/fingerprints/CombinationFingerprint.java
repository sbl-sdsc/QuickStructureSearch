package org.rcsb.fingerprints;

import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Point3d;

/**
 * This class generates a fingerprint (signature) for protein chains for null hypothesis testing.
 * The null hypothesis H0 is: all protein fragments of equal length are considered equal for 
 * distinguishing similar from dissimilar proteins. To reject the null hypothesis, the statistical 
 * significance (i.e. 95% or 99% likelihood) of an alternative fingerprint method must be demonstrated.
 * See for example: 
 * https://explorable.com/null-hypothesis
 * http://www.null-hypothesis.co.uk/science/item/what_is_a_null_hypothesis/
 *
 * @author Peter Rose
 */
public class CombinationFingerprint implements GenericFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private GenericFingerprint fingerprint1;
	private GenericFingerprint fingerprint2;
	private double weight1;
	private double weight2;	
  
    /**
     * Default constructor uses default parameters
     */
    public CombinationFingerprint(GenericFingerprint fingerprint1, GenericFingerprint fingerprint2, double weight1, double weight2) {
    	this.fingerprint1 = fingerprint1;
    	this.fingerprint2 = fingerprint2;
    	this.weight1 = weight1;
    	this.weight2 = weight2;
    }
    
    public String getName() {
    	return fingerprint1.getName()  + fingerprint2.getName();
    }
    
    /**
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
    public double[] getFingerprint(Point3d[] points) {
       double[] f1 = fingerprint1.getFingerprint(points);
       double[] f2 = fingerprint2.getFingerprint(points);
 //      System.out.println("f1: " + Arrays.toString(f1));
 //      System.out.println("f2: " + Arrays.toString(f2));
       double sum1 = sum(f1);
       double sum2 = sum(f2);
       
       double relWeight = 1;
       if (sum1 > 0 && sum2 > 0) {
           relWeight = sum2/sum1;
       }
//       System.out.println("sums: " + sum1 + ", " + sum2);
       f1 = scale(f1, weight1*relWeight);
       f2 = scale(f2, weight2);
 //      System.out.println("wf1: " + Arrays.toString(f1));
 //      System.out.println("wf2: " + Arrays.toString(f2));

       double[] f = Arrays.copyOf(f1, f1.length+f2.length);
       System.arraycopy(f2, 0, f, f1.length, f2.length);
//       System.out.println("f: " + Arrays.toString(f));
       return f;
	}
    
    private double sum(double[] values) {
    	double sum = 0;
    	for (double v: values) {
    		sum += v;
    	}
        return sum;
    }
    
    private double[] scale(double[] values, double scaleFactor) {
    	for (int i = 0; i < values.length; i++) {
    		values[i] = values[i] * scaleFactor;
    	}
    	return values;
    }
}
