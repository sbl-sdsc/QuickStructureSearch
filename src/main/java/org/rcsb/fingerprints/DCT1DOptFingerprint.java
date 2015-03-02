package org.rcsb.fingerprints;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.vecmath.Point3d;

/**
 * Creates a fingerprint for protein fragments that represents the overall protein chain structure
 * based on the one dimensional discrete cosine transform 
 * 
 * @author Alan Yeung, Peter Rose
 */
public class DCT1DOptFingerprint implements GenericFingerprint, LinearFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private DiscreteCosineTransform transform = new DiscreteCosineTransform();
	// parameter set 1 (9,1,1,9,15,5,5,80 sens=0.976, spec=0.561, positive hashcodes only)
	// parameter set 2 (9,1,2,9,20,4,4,80 sens=0.991, spec=0.501, posiitve hashcodes only)
	// parameter set 3 (8,1,4,5,13,1,1,80 sens=0.963, spec=0.551, positive and negative hashcode (abs value)), better w/ dim=200: 0.94/0.64
	// parameter set 4 (9,1,1,5,22,1,1,80 sens=0.987, spec=0.488, --------------- " ------------------------), better w/ dim=200: 0.97/0.60
	// parameter set 5 (10,1,1,6,5,1,1,80  sens=0.971, spec=0.544, --------------- " ------------------------)
	// parameter set 6 (8,1,1,8,9,1,1,200 sens=0.966, spec=0.552, ")
	// parameter set 7 (9,1,1,5,25,2,1,200 sens=0.990, spec=0.418, ")
	private int length = 10;
	private int skip = 1;
	private int minDistance = 1;
	private int terms = 6;
	private double dcCoefficient = 5; 
	private double acCoefficientHigh = 1;
	private double acCoefficientLow = 1;
	private double lNorm = 2.0;
	private int dimensions = 200;

    private int[][] distancePairs;
    private double[] coefficients;
	
	public static void main(String args[]) {
		DCT1DOptFingerprint f = new DCT1DOptFingerprint();
		f.getDistancePairs();
	}
	public DCT1DOptFingerprint(){
		distancePairs = getDistancePairs();
		coefficients = getQuantizationCoefficients();
	}

// parameter sets: 1,2
//	public DCT1DOptFingerprint(int randomSeed){
//		Random r = new Random(randomSeed);
//		this.length = 7 + r.nextInt(10);
//		this.minDistance = 1 + r.nextInt(4);
//		this.terms = 5 + r.nextInt(10);
//		if (this.terms < this.length) {
//			this.terms = this.length;
//		}
//		this.dcCoefficient = 10 + r.nextInt(21);
//		this.acCoefficientHigh = 3 + r.nextInt(11);
//		this.acCoefficientLow = 2 + r.nextInt(11);
//		if (acCoefficientHigh < acCoefficientLow) {
//			double tmp = acCoefficientLow;
//			acCoefficientLow = acCoefficientHigh;
//			acCoefficientHigh = tmp;
//		}
//	
//		distancePairs = getDistancePairs();
//		coefficients = getQuantizationCoefficients();
//	}

// parameter sets: 3,4,5
//	public DCT1DOptFingerprint(int randomSeed){
//		Random r = new Random(randomSeed);
//		this.length = 7 + r.nextInt(4);
//		this.minDistance = 1 + r.nextInt(4);
//		this.terms = 5 + r.nextInt(4);
//		this.dcCoefficient = 5 + r.nextInt(21);
//		this.acCoefficientHigh = 1 + r.nextInt(11);
//		this.acCoefficientLow = 1 + r.nextInt(11);
//		if (acCoefficientHigh < acCoefficientLow) {
//			double tmp = acCoefficientLow;
//			acCoefficientLow = acCoefficientHigh;
//			acCoefficientHigh = tmp;
//		}
//	
//		distancePairs = getDistancePairs();
//		coefficients = getQuantizationCoefficients();
//	}
	// parameter sets: 6, 7
	public DCT1DOptFingerprint(int randomSeed){
		Random r = new Random(randomSeed);
		this.length = 8 + r.nextInt(3);
		this.minDistance = 1;
		this.terms = 3 + r.nextInt(6);
		this.dcCoefficient = 5 + r.nextInt(21);
		this.acCoefficientHigh = 1 + r.nextInt(11);
		this.acCoefficientLow = 1 + r.nextInt(11);
		if (acCoefficientHigh < acCoefficientLow) {
			double tmp = acCoefficientLow;
			acCoefficientLow = acCoefficientHigh;
			acCoefficientHigh = tmp;
		}
		this.dimensions = 200;
	
		distancePairs = getDistancePairs();
		coefficients = getQuantizationCoefficients();
	}
	
	public float[] getParameters() {
		float[] parameters = {
				this.length,
				this.skip,
				this.minDistance,
				this.terms,
				(float) this.dcCoefficient,
				(float) acCoefficientHigh,
				(float) acCoefficientLow,
				this.dimensions
				};
		return parameters;
	}
	
	public DCT1DOptFingerprint(int length, int skip, int minDistance, int terms, double dcCoefficient, double acCoefficientHigh, double acCoefficientLow, int dimensions){
		this.length = length;
		this.skip = skip;
		this.minDistance = minDistance;
		this.terms = terms;
		this.dcCoefficient = dcCoefficient;
		this.acCoefficientHigh = acCoefficientHigh;
		this.acCoefficientLow = acCoefficientLow;
		this.dimensions = dimensions;

		distancePairs = getDistancePairs();
		coefficients = getQuantizationCoefficients();
	}
	
	public double[] getFingerprint(Point3d[] coords) {
		double[] features = new double[this.dimensions];
		
		if(coords.length-this.length < 0){
    		return features;
    	}

		double[] distances = new double[distancePairs.length];	
    	
    	for(int i = 0; i < coords.length-this.length+1; i+= skip){
    		if (hasGaps(coords, i)){
    			continue; //Skip if gap exists
    		}
    		
    		distances = getDistanceMatrix(coords, i);
    		double[] dct = transform.getTransform(distances);
 //   		System.out.println("dct: " + Arrays.toString(dct));
            int[] dctCoeff = DiscreteCosineTransform.quantize(dct, coefficients);
 //           double[] recovered = DiscreteCosineTransform.dequantize(dctCoeff, coefficients, distances.length);
//            double[] idct = transform.getInverseTransform(recovered);
//            System.out.println(Arrays.toString(dctCoeff));
//            System.out.println("d:  " + Arrays.toString(distances));
//            System.out.println("r:  " + Arrays.toString(idct));
//            System.out.println("rmsd: " + DiscreteCosineTransform.rmsd(distances, idct));
 //          System.out.println("quantized: " + Arrays.toString(values));
    		
    		// calculate unique hash code for fragment
    		int hashCode = Arrays.hashCode(dctCoeff);
    		
    		// map hashCode into feature vector of fixed size by hashing.
    		// This is called the "hashing trick". See http://en.wikipedia.org/wiki/Feature_hashing
    		features[Math.abs(hashCode) % this.dimensions]++;
    	}

//    	Set<Integer> set = new HashSet<Integer>();
//    	for (int h: getLinearFingerprint(coords)) {
//    		set.add(h);
//    	}
 //   	System.out.println("Set: " + set.size());
//    	System.out.println(Arrays.toString(getLinearFingerprint(coords)));
    	return features;
	}
	
	public int[] getLinearFingerprint(Point3d[] coords) {
		int[] fingerprint = new int[coords.length-this.length+1];
		
		if(coords.length-this.length < 0){
    		return fingerprint;
    	}

		double[] distances = new double[distancePairs.length];	
    		
    	for(int i = 0; i < coords.length-this.length+1; i+= skip){
    		if (hasGaps(coords, i)){
    			continue; //Skip if gap exists
    		}
    		
    		distances = getDistanceMatrix(coords, i);
    		double[] dct = transform.getTransform(distances);
 //   		System.out.println("dct: " + Arrays.toString(dct));
            int[] dctCoeff = DiscreteCosineTransform.quantize(dct, coefficients);
 //           double[] recovered = DiscreteCosineTransform.dequantize(dctCoeff, coefficients, distances.length);
//            double[] idct = transform.getInverseTransform(recovered);
//            System.out.println(Arrays.toString(dctCoeff));
//            System.out.println("d:  " + Arrays.toString(distances));
//            System.out.println("r:  " + Arrays.toString(idct));
//            System.out.println("rmsd: " + DiscreteCosineTransform.rmsd(distances, idct));
 //          System.out.println("quantized: " + Arrays.toString(values));
    		
    		// calculate unique hash code for fragment
    		int hashCode = Arrays.hashCode(dctCoeff);
  //  		if (hashCode < 0) {
    //			continue;
  //  		}
    		
    		// map hashCode into feature vector of fixed size by hashing.
    		// This is called the "hashing trick". See http://en.wikipedia.org/wiki/Feature_hashing
 //   		fingerprint[i] = hashCode % this.dimensions;
    		fingerprint[i] = Math.abs(hashCode);
    	}

    	return fingerprint;
	}
	
	/**
	 * Returns a string representing the name of this fingerprint
	 * @return name of this fingerprint and the length of fragment used
	 */
	public String getName(){
		return this.getClass().getSimpleName() + "_L" + this.length + "T" + this.terms;  
	}
	
	/**
	 * 
	 * @param coords
	 * @param index
	 * @return
	 */
	private double[] getDistanceMatrix(Point3d[] coords, int index){
		double[] distance = new double[distancePairs.length];

		for (int i = 0; i < distancePairs.length; i++){
			int j = distancePairs[i][0] + index;
			int k = distancePairs[i][1] + index;
			distance[i] = coords[j].distance(coords[k]);
//			double dl = Math.pow(Math.abs(coords[j].x-coords[k].x), lNorm) + 
//					Math.pow(Math.abs(coords[j].y-coords[k].y), lNorm) +
//					Math.pow(Math.abs(coords[j].z-coords[k].z), lNorm);
//			distance[i] = Math.abs(Math.pow(dl, 1.0/lNorm));
		}
		return distance;
	}
	
	/**
	 * Returns true if there is a gap between the C alpha atoms
	 * within a fragment. Note, the first position is not checked, 
	 * since the gap info is for a gap between residue i and i + 1.
	 * @param gaps true if there is a gap in the fragment
	 * @param index start residue of fragment
	 * @return
	 */
	private boolean hasGaps(Point3d[] coords, int index) {
		for (int i = index; i < index+this.length; i++) {
			if (coords[i] == null) {
				return true;
			}
		}
		return false;
	}
	
	private int[][] getDistancePairs() {
		int n = 0;
		for (int d = this.length-1; d >= minDistance; d--) {
			for (int i = 0; i < this.length-d; i++) {
				n++;
			}
		}
		int[][] pairs = new int[n][2];
		
		int j = 0;
		for (int d = this.length-1; d >= minDistance; d--) {
			for (int i = 0; i < this.length-d; i++) {
				pairs[j][0] = i;
				pairs[j][1] = i + d;
				j++;
			}
		}
		return pairs;
	}
	
	private double[] getQuantizationCoefficients() {
		double[] coefficients = new double[this.terms];
		coefficients[0] = dcCoefficient;
		
		System.out.println("terms:  "+ this.terms);
		double alpha = 1.0/(this.terms-2);
		for (int i = 0; i < this.terms-1; i++) {
			coefficients[i+1] = acCoefficientHigh * (1-alpha*i) + acCoefficientLow * alpha*i;
			System.out.println("alpha: " + alpha*i);
		}
		return coefficients;
	}
}
