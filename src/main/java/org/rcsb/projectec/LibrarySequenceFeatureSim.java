package org.rcsb.projectec;

import java.io.Serializable;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.commons.lang.ArrayUtils;
import org.rcsb.project3.SequenceFeatureInterface;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for storing the list of indicies resulting from LibraryFingerprint
 * No similarity for this, may be used for Jaccard Index
 * 
 * @author Emilia Copic
 */
public class LibrarySequenceFeatureSim implements SequenceFeatureInterface<Integer>, Serializable {

	private static final long serialVersionUID = 1L;
	List<Point3d[]> library;
	double[][] rmsdArray;
	private int[] libSequence;
	private Point3d[] coords;
	// Some setting for calculate similarity score 
	double rmsdThreshold;
	private double gap = -1;
    /**
     * Constructor that will store an int array as EndToEndDistance feature
     * @param EndToEndSequence
     */
	public LibrarySequenceFeatureSim(int[] LibrarySequence, double[][] rmsdArray, double rmsdThreshold) {
		this.libSequence = LibrarySequence;
		this.rmsdThreshold = rmsdThreshold;
		this.rmsdArray = rmsdArray;
		
	}
	
	public LibrarySequenceFeatureSim(int[] EndToEndSequence, Point3d[] coords,  List<Point3d[]> lib) {
		this.libSequence = EndToEndSequence;
		this.coords = coords;
		this.library = lib;
	}
		
	@Override
	public double similarity(SequenceFeatureInterface<Integer> sequence2, int i, int j) {
		// check NaN as gap

		if (this.get(i) == null || sequence2.get(j) == null){
			return gap;
			//I might have gotten rid of all nulls...
		}
		
		else{
			
		double rmsd = rmsdArray[this.get(i)][sequence2.get(j)];
		double score;
		score = 1 - 2*(rmsd/rmsdThreshold);
//		if(rmsd>2*rmsdThreshold)
//		score = -2;
//		
//		else if(rmsd > 1.5*rmsdThreshold)
//			score = -1;
//		
//		else if(rmsd >rmsdThreshold)
//			score = -.5;
//		
//		else if(rmsd > .7*rmsdThreshold)
//			score = 1;
//		else 
//			score = 2;
			//(rmsd < rmsdThreshold)
		
		return score;
		}
	}

	@Override
	public boolean identity(SequenceFeatureInterface<Integer> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return false;
		}
		// check identity
		else 
			return this.get(i).equals(sequence2.get(j));
	}

	@Override
	public Integer get(int index) {
		if (Double.isNaN((double)libSequence[index]))
			return null;
		else
			return libSequence[index];
	}

	@Override
	public int length() {
		return libSequence.length;
	}

	@Override
	public Integer[] getSequence() {
		return ArrayUtils.toObject(libSequence);
	}

	@Override
	public String toString(int index) {
		if (this.get(index) == null)
			return "NULL";
		else 			
			return this.get(index).toString();
	}

	@Override
	public double todouble(int index) {
		return (double)(int)this.get(index);
	}

	@Override
	public Point3d[] getCoords() {
		return coords;
	}
}
