package org.rcsb.fingerprints;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;

import javax.vecmath.Point3d;

/**
 * Creates a fingerprint for protein fragments that represents the overall protein chain structure
 * based on the one dimensional discrete cosine transform 
 * 
 * @author Alan Yeung, Peter Rose
 */
public class DCT1DLinearFingerprint implements LinearFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private DiscreteCosineTransform transform = new DiscreteCosineTransform();
	// parameter set 1 (9,1,1,9,15,5,5,80 sens=0.976, spec=0.561, positive hashcodes only)
	// parameter set 2 (9,1,2,9,20,4,4,80 sens=0.991, spec=0.501, posiitve hashcodes only)
	// parameter set 3 (8,1,4,5,13,1,1,80 sens=0.963, spec=0.551, positive and negative hashcode (abs value)), better w/ dim=200: 0.94/0.64
	// parameter set 4 (9,1,1,5,22,1,1,80 sens=0.987, spec=0.488, --------------- " ------------------------), better w/ dim=200: 0.97/0.60
	// parameter set 5 (10,1,1,6,5,1,1,80  sens=0.971, spec=0.544, --------------- " ------------------------)
	// parameter set 6 (8,1,1,8,9,1,1,200 sens=0.966, spec=0.552, ")
	// parameter set 7 (9,1,1,5,25,2,1,200 sens=0.990, spec=0.418, ")
	// parameter set 8 (10,1,1,8,8,9,5,200 sens=0.982, sens=0.541, with Levenshtein score)
	// parameter set 9 (10,1,1,3,12,4,1,200 sens=0.926, spec=0.930, opt. w/ Levenshtein score)
	// parameter set10 (10,1,1,3,17,4,2,200 sens=0.946, spec=0.903, ------------ " -----------
	// parameter set11 (9,1,1,3,18,6,1,200  sens=0.910, spec=0.944, ------------ " -----------
	// parameter set12 (9,1,1,3,10.70,3.64,1.34,200 sens=0.975, spec=0.833, ----- " ----------
	// parameter set13 (9,1,1,3,10.93,3.35,1.05,200 sens=0.932, spec=0.935
	
	private int length = 9;
	private int skip = 1;
	private int minDistance = 1;
	private int terms = 3;
	private double dcCoefficient = 10.70; 
	private double acCoefficientHigh = 3.64;
	private double acCoefficientLow = 1.05;
	private double pNorm = 2.0;
	private int dimensions = 200;

    private int[][] distancePairs;
    private double[] coefficients;
	
	public static void main(String args[]) {
		DCT1DLinearFingerprint f = new DCT1DLinearFingerprint();
		f.getDistancePairs();
	}
	public DCT1DLinearFingerprint(){
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
	public DCT1DLinearFingerprint(int randomSeed){
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
	
	public DCT1DLinearFingerprint(int length, int skip, int minDistance, int terms, double dcCoefficient, double acCoefficientHigh, double acCoefficientLow, int dimensions){
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
    		
    		distances = getDistancePairs(coords, i, distances);
 //   		System.out.println("d: " + Arrays.toString(distances));
    		double[] dct = transform.getTransform(distances);
 //   		System.out.println("dct: " + Arrays.toString(dct));
            int[] dctCoeff = DiscreteCosineTransform.quantize(dct, coefficients);
 //           double[] recovered = DiscreteCosineTransform.dequantize(dctCoeff, coefficients, distances.length);
 //           double[] idct = transform.getInverseTransform(recovered);
//            System.out.println(Arrays.toString(dctCoeff));
 //           System.out.println("d:  " + Arrays.toString(distances));
 //           System.out.println("r:  " + Arrays.toString(idct));
 //           System.out.println("rmsd: " + DiscreteCosineTransform.rmsd(distances, idct));
 //          System.out.println("quantized: " + Arrays.toString(values));
    		
    		// calculate unique hash code for fragment
    		int hashCode = Arrays.hashCode(dctCoeff);
    		
    		// map hashCode into feature vector of fixed size by hashing.
    		// This is called the "hashing trick". See http://en.wikipedia.org/wiki/Feature_hashing
 //   		fingerprint[i] = hashCode % this.dimensions;
 //   		fingerprint[i] = Math.abs(hashCode);
    		fingerprint[i] = hashCode;
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
