package org.rcsb.projectm;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;

import javax.vecmath.Point3d;

import org.rcsb.fingerprints.DiscreteCosineTransform;
import org.rcsb.project3.DCT1DSequenceFeature;
import org.rcsb.project3.SequenceFingerprint;

/**
 * Creates a fingerprint for protein fragments that represents the overall protein chain structure
 * based on the one dimensional discrete cosine transform 
 * 
 * @author Alan Yeung, Peter Rose
 */
public class DCT1DSequenceFingerprint implements SequenceFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private DiscreteCosineTransform transform = new DiscreteCosineTransform();
	
	private int length = 9;
	private int skip = 1;
	private int minDistance = 1;
	private int terms = 3;
	private double dcCoefficient = 10.70; 
	private double acCoefficientHigh = 3.64;
	private double acCoefficientLow = 1.05;
	private double pNorm = 2.0;
	private int dimensions = 200;
	
	// Flag for if this object is created with settings
	private boolean settingFlag = false;
	private double gap, match, mismatch;

    private int[][] distancePairs;
    private double[] coefficients;
	
	public static void main(String args[]) {
		DCT1DSequenceFingerprint f = new DCT1DSequenceFingerprint();
		f.getDistancePairs();
	}
	
	public DCT1DSequenceFingerprint(){
		settingFlag = false;
		distancePairs = getDistancePairs();
		coefficients = getQuantizationCoefficients();
	}
	
	public DCT1DSequenceFingerprint(double gap, double match, double mismatch){
		distancePairs = getDistancePairs();
		coefficients = getQuantizationCoefficients();
		settingFlag = true;
		this.gap = gap;
		this.match = match;
		this.mismatch = mismatch;
	}

	public DCT1DSequenceFingerprint(int randomSeed){
		settingFlag = false;
		Random r = new Random(randomSeed);
		this.length = 8 + r.nextInt(9);
		this.skip = 1 + r.nextInt(8);
		this.minDistance = 1 + r.nextInt(2);
		this.terms = 3 + r.nextInt(3);
		this.dcCoefficient = 10 + r.nextDouble()*10;
		this.acCoefficientHigh = 1 + r.nextDouble()*5;
		this.acCoefficientLow = r.nextDouble()*3;
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
	
	public DCT1DSequenceFingerprint(int length, int skip, int minDistance, int terms, double dcCoefficient, double acCoefficientHigh, double acCoefficientLow, int dimensions){
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
	
	public DCT1DSequenceFeature getFingerprint(Point3d[] coords) {
		int[] fingerprint = new int[coords.length-this.length+1];
		
		if(coords.length-this.length < 0){
			if (settingFlag)
				return new DCT1DSequenceFeature(fingerprint,gap,match,mismatch);
			else
				return new DCT1DSequenceFeature(fingerprint);
    	}

		double[] distances = new double[distancePairs.length];	
    		
    	for(int i = 0; i < coords.length-this.length+1; i+= skip){
    		if (hasGaps(coords, i)){
    			continue; //Skip if gap exists
    		}
    		
    		distances = getDistancePairs(coords, i, distances);
    		double[] dct = transform.getTransform(distances);
            int[] dctCoeff = DiscreteCosineTransform.quantize(dct, coefficients);
    		
    		// calculate unique hash code for fragment
    		int hashCode = Arrays.hashCode(dctCoeff);
    		fingerprint[i] = hashCode;
    	}
    	
    	if (settingFlag)
			return new DCT1DSequenceFeature(fingerprint,gap,match,mismatch);
		else
			return new DCT1DSequenceFeature(fingerprint);
	}
	
	/**
	 * Returns a string representing the name of this fingerprint
	 * @return name of this fingerprint and the length of fragment used
	 */
	public String getName(){
		return this.getClass().getSimpleName() + "_L" + this.length + "T" + this.terms;  
	}
	
	/**
	 * Returns distances between pairs of C-alpha atoms
	 * in a fragment
	 * @param coords
	 * @param index
	 * @return
	 */
	private double[] getDistancePairs(Point3d[] coords, int index, double[] distances) {
		for (int i = 0; i < distancePairs.length; i++){
			int j = distancePairs[i][0] + index;
			int k = distancePairs[i][1] + index;
			distances[i] = coords[j].distance(coords[k]);
		}
		return distances;
	}
	/**
	 * Returns the L^p Norm for C-alpha atom pairs in 
	 * a fragment
	 * @param coords
	 * @param index
	 * @return
	 */
	@SuppressWarnings("unused")
	private double[] getLpNormPairs(Point3d[] coords, int index, double[] distances) {
		for (int i = 0; i < distancePairs.length; i++){
			int j = distancePairs[i][0] + index;
			int k = distancePairs[i][1] + index;
			double dl = Math.pow(Math.abs(coords[j].x-coords[k].x), pNorm) + 
					Math.pow(Math.abs(coords[j].y-coords[k].y), pNorm) +
					Math.pow(Math.abs(coords[j].z-coords[k].z), pNorm);
			distances[i] = Math.abs(Math.pow(dl, 1.0/pNorm));
		}
		return distances;
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
	
	/**
	 * Returns pairs of indices that specify distance pairs within
	 * a fragment. The pairs are sorted from longest to shortest
	 * distances
	 */
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
	
	/**
	 * Returns a list of quantization coefficients to convert
	 * DCT coefficients to integer values
	 */
	private double[] getQuantizationCoefficients() {
		double[] coefficients = new double[this.terms];
		
		// the first term is the most important term
		coefficients[0] = dcCoefficient;
		
		// for the remaining terms the values are scaled down by a ramp
		// starting at the acCoefficientHigh and acCoefficientLow
		double alpha = 1.0/(this.terms-2);
		for (int i = 0; i < this.terms-1; i++) {
			coefficients[i+1] = acCoefficientHigh * (1-alpha*i) + acCoefficientLow * alpha*i;
		}
		return coefficients;
	}
}
