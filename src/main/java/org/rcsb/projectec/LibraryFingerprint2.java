package org.rcsb.projectec;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.rcsb.project3.EndToEndDistanceSequenceFeature;
import org.rcsb.project3.SequenceFeatureInterface;
import org.rcsb.project3.SequenceFingerprint;
import org.rcsb.structuralAlignment.SuperPosition;
import org.rcsb.structuralAlignment.SuperPositionQCP;

/**
 * Takes a library and identifies a protein as a list of fragment indices
 *  from the library.
 * @author Emilia Copic
 */

public class LibraryFingerprint2 implements SequenceFingerprint, Serializable {
	// will take a protein, convert it to indices of library
	private int length = 8;
	List<Point3d[]> library;
	double[][] rmsdArray;
	private SuperPositionQCP qcp = new SuperPositionQCP(true);
	private double rmsdThreshold;

	/**
	 * Constructor with all parameters
	 * 
	 * @param lib
	 *            library with all of the archetype fragments
	 * @param rmsdArray 
	 * @param rmsdThreshold
	 *            threshold it will check the library against - make this the
	 *            same threshold as the library
	 * 
	 */
	public LibraryFingerprint2 (List<Point3d[]> lib, double[][] rmsdArray, double rmsdThreshold) {
		this.library = lib;
		this.rmsdThreshold = rmsdThreshold;
		length = library.get(0).length;
		this.rmsdArray = rmsdArray;
	}

	public int getLength() {
		return length;
	}

	public String getName() {
		return this.getClass().getSimpleName() + "_L" + this.length + "B";
	}

	/**
	 * Returns a fingerprint for the given chain.
	 * 
	 * @param coords
	 *            coordinates of a macromolecule fragment
	 * @return fingerprint
	 */
	@Override
	public LibrarySequenceFeatureSim getFingerprint(Point3d[] coords){
		// split protein into fragments
		List<Point3d[]> fragments = new ArrayList<>();
		for (int i = 0; i < coords.length - length + 1; i++) {
		Point3d[] fragment = Arrays.copyOfRange(coords, i, i + length);
			boolean result = true;
			for (int j = 0; j < fragment.length; j++) {
				result = result && (fragment[j] != null);
			}
			//if no nulls, add
			if (result) {	
			fragments.add(fragment);
		}
		}

		// compare fragments to each one in library
		//then get list of indices the fragments correspond to
		List<Integer> typeIndices = new ArrayList<>();
		for (int j = 0; j < fragments.size(); j++) {
			
			Point3d[] cFragment = SuperPosition.clonePoint3dArray(fragments.get(j));
			SuperPositionQCP.center(cFragment);
			// compare each archetype with the new fragment
			// for (Point3d[] archetype: library) {

			for (int i = 0; i < library.size(); i++) {
				qcp.set(library.get(i), cFragment);
				double rmsd = qcp.getRmsd();

				if (rmsd < rmsdThreshold) {
					typeIndices.add(i);
					break;
					//not sure what to do if for some reason they don't match...
				}
			}

		}
		//make typeIndices into an int array
		int[] typeIndicesArray = new int[typeIndices.size()];
		for (int i = 0; i < typeIndices.size(); i++) {
			typeIndicesArray[i] = typeIndices.get(i).intValue();
		}
		
		return new LibrarySequenceFeatureSim(typeIndicesArray, rmsdArray, rmsdThreshold);
		
		
	

	}

}
