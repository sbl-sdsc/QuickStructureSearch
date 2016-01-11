package org.rcsb.fingerprints;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

/**
 * This class generates a fingerprint (signature) for protein chains based quartets
 * of points that span a tetrahedron in a protein chains
 *
 * @author Peter Rose
 */
public class TetrahedronFingerprint implements GenericFingerprint, Serializable {
	private static final long serialVersionUID = 1L;
	private int length = 9;
	private double binSize = 3.0;
	private int dimensions = 80; // 40: 821056, 60: 469370, 80: 420455
	private int increment = 1;
//	private List<Integer> list = Arrays.asList(4,8,16,32,64);
	private List<Integer> list = Arrays.asList(8,16,32,64);
//	private List<Integer> list = Arrays.asList(8,16,24,32,40,48,56,64);

	/**
	 * Default constructor uses default parameters
	 */
	public TetrahedronFingerprint() {}

	/**
	 * Constructor with all parameters
	 * @param length fragment length
	 */
	public TetrahedronFingerprint (int length, double binSize) {
		this.length = length;
		this.binSize = binSize;
	}

	public String getName() {
		return this.getClass().getSimpleName() + "_L" + this.length + "B" + this.binSize;
	}

	/**
	 * Returns a fingerprint for the given chain. 
	 * @param coords coordinates of a macromolecule fragment
	 * @return fingerprint
	 */
	public double[] getFingerprint(Point3d[] points) {
		double scale = 1/binSize;

		int maxIndex = 0;
		int[][] dm = new int[points.length][points.length];

		for (int i = 0; i < points.length-length-1; i+=increment) {
			for (int j = i + length; j < points.length; j+=increment) {
				int index = 0;
				if (points[i] != null && points[j] != null) {
					index = (int)Math.floor(scale*points[i].distance(points[j]));
					if (!list.contains(index)) {
						index = 0;
					}
					//				System.out.println(index + ": " + points[i].distance(points[j]));
					maxIndex = Math.max(maxIndex,  index);
				}
				dm[i][j] = index;
				dm[j][i] = index;		
			}
		}

		double features[] = new double[dimensions];

		int[] distances = new int[6];
		for (int i = 0; i < points.length-3*length-3; i+=increment) {
			for (int j = i + length; j < points.length-2*length-2; j+=increment) {
				if (dm[i][j] != 0) {
					for (int k = j + length; k < points.length-length-1; k+=increment) {
						if (dm[i][k] != 0 && dm[j][k] != 0 && dm[i][j] == dm[i][k]) {
							for (int l = k + length; l < points.length; l+=increment) {
//								if (dm[i][l] != 0 && dm[i][j] == dm[i][l] && dm[j][l] != 0 && dm[k][l] != 0) {
									if (dm[i][l] != 0 && dm[j][l] != 0 && dm[k][l] != 0) {
									distances[0] = dm[i][j];
									distances[1] = dm[i][k];
									distances[2] = dm[i][l];
									distances[3] = dm[j][k];
									distances[4] = dm[j][l];
									distances[5] = dm[k][l];
					//				Arrays.sort(distances);
									// calculate unique hash code for fragment
									int hashCode = Arrays.hashCode(distances);

									// map hashCode into feature vector of fixed size by hashing.
									// This is called the "hashing trick". See http://en.wikipedia.org/wiki/Feature_hashing
									int index = hashCode % this.dimensions;
									if (index >= 0) {
										features[index]++;
									}
								}
							}
						}
					}
				}
			}
		}
	//	System.out.println(maxIndex);
	//		System.out.println(Arrays.toString(features));

		return features;
	}
}
