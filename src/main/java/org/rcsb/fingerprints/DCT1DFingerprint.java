package org.rcsb.fingerprints;

import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Point3d;

/**
 * Creates a fingerprint for protein fragments that represents the overall protein chain structure
 * based on the one dimensional discrete cosine transform 
 * 
 * @author Alan Yeung, Peter Rose
 */
public class DCT1DFingerprint implements GenericFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private int length = 8;
	private int terms = 6;
	private int binWidth = 15;
	private int dimensions = 40;

//	private static final int[][] distancePairs = 
//	{{0,7},{0,6},{1,7},{0,5},{1,6},{2,7},{0,4},{1,5},{2,6},{3,7},{0,3},{1,4},{2,5},{3,6},{4,7},{0,2},{1,3},{2,4},{3,5},{4,6},{5,7}};
	private int[][] distancePairs = getDistancePairs();
	private double[][] cosMatrix = initializeCosMatrix(distancePairs.length);
	
	
	public static void main(String args[]) {
		DCT1DFingerprint f = new DCT1DFingerprint();
		f.getDistancePairs();
	}
	public DCT1DFingerprint(){
	}
	
	public DCT1DFingerprint(int dimensions){
		this.dimensions = dimensions;
	}
	
	public DCT1DFingerprint(int length, int dimensions){
		this.length = length;
		this.dimensions = dimensions;
	}
	
	public double[] getFingerprint(Point3d[] coords) {
		double[] features = new double[this.dimensions];
		
		if(coords.length-this.length-1 <= 0){
    		return features;
    	}

		double[] distances = new double[distancePairs.length];	
    	
    	for(int i = 0; i < coords.length-this.length-1; i++){
    		if (hasGaps(coords, i)){
    			continue; //Skip if gap exists
    		}
    		
    		distances = getDistanceMatrix(coords, i);
    		double[] dct = getDCT(distances); 
            int[] values = quantize(dct);
 //          System.out.println("quantized: " + Arrays.toString(values));
    		
    		// calculate unique hash code for fragment
    		int hashCode = Arrays.hashCode(values);
    		if (hashCode < 0) {
    			continue;
    		}
    		
    		// map hashCode into feature vector of fixed size by hashing.
    		// This is called the "hashing trick". See http://en.wikipedia.org/wiki/Feature_hashing
    		features[hashCode % this.dimensions]++;
    	}

    	return features;
	}
	
	/**
	 * Returns a string representing the name of this fingerprint
	 * @return name of this fingerprint and the length of fragment used
	 */
	public String getName(){
		return this.getClass().getSimpleName() + "_L" + this.length + "T" + this.terms + "B" + this.binWidth;  
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
		}
		return distance;
	}
    
    /**
     * Takes in a 2D array and returns the discrete cosine transform of it
     */
    private double[] getDCT (double[] distances){
    	double[] dct = new double[distances.length];
    	
    	for (int u = 0; u < distances.length; u++){
    		   dct[u] = getSummationValue(distances, u) * getAlphaConstant(u);	
    	}	
    	return dct;
    }
	
	/**
	 * Returns alpha constant used in caluclation of discrete cosine transform
	 * @param index
	 * @return
	 */
	private static double getAlphaConstant(int index){
		if(index == 0){
			return 1/Math.sqrt(2);
		}
		else
			return 1;
	}
	
	/**
	 * Returns the summation value part of the discrete cosine transform formula
	 * @param distanceMatrix the matrix to get values from for the summation
	 * @param u horizontal spatial frequency of dct transformed matrix
	 * @param v vertical spatial frequency of dct transformed matrix
	 * @return a double with the value of the summation part of the dct formula
	 */
	private double getSummationValue(double[] distances, int u) {
		double sum = 0;
		
		for (int x = 0; x < distances.length; x++){
			sum += distances[x] * cosMatrix[x][u];
		}
		
		return sum;
	}
	
	/**
	 * Returns the values for the cos calculations in the DCT summation part.
	 * @return the cos calculations with the x/y value on the horizontal,
	 *    and the u/v values on the vertical
	 */
	private static double[][] initializeCosMatrix(int terms) {
		double[][] cosMatrix= new double[terms][terms];
		
		for(int xy = 0; xy < terms; xy++){
			for(int uv = 0; uv < terms; uv++){
				cosMatrix[xy][uv] = Math.cos((( 2*xy + 1) * uv*Math.PI)/ (2 * terms));
			}
		}	
		return cosMatrix;	
	}
		
	private int[] quantize(double[] dct) {
		int[] values = new int[this.terms];

		// quantize the DC component
		values[0] = (int)Math.round(dct[0]/(3*this.binWidth));
		
		// quantize the AC components
		for (int i = 1; i < values.length; i++) {
			values[i] = (int) Math.round(dct[i]/this.binWidth);
		}
		
		return values;
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
		for (int i = 0; i < this.length-1; i++) {
			n += i;
		}
		int[][] pairs = new int[n][2];
		
		int j = 0;
		for (int d = this.length-1; d > 1; d--) {
			for (int i = 0; i < this.length-d; i++) {
				pairs[j][0] = i;
				pairs[j][1] = i + d;
				j++;
			}
		}
		return pairs;
	}
}
