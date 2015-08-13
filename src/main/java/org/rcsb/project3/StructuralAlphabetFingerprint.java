package org.rcsb.project3;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
/**
 * This class transform protein's psi and phi angle into unique blocks.
 * In order to use this fingerPrint, may need to use sequence file that contain CA, N and C atoms
 * 
 * @author Chris Li
 */
public class StructuralAlphabetFingerprint implements SequenceFingerprint, Serializable {

	private static final long serialVersionUID = 1L;
	// most frequent blocks
	double[][] references = {
			{ 41.14,   75.53,  13.92,  -99.80,  131.88,  -96.27, 122.08,  -99.68},
			{108.24,  -90.12, 119.54,  -92.21,  -18.06, -128.93, 147.04,  -99.90},
			{-11.61, -105.66,  94.81, -106.09,  133.56, -106.93, 135.97, -100.63},
			{141.98, -112.79, 132.20, -114.79,  140.11, -111.05, 139.54, -103.16},
			{133.25, -112.37, 137.64, -108.13,  133.00,  -87.30, 120.54,   77.40},
			{116.40, -105.53, 129.32,  -96.68,  140.72,  -74.19, -26.65,  -94.51},
			{  0.40,  -81.83,   4.91, -100.59,   85.50,  -71.65, 130.78,   84.98},
			{119.14, -102.58, 130.83,  -67.91,  121.55,   76.25,  -2.95,  -90.88},
			{130.68,  -56.92, 119.26,   77.85,   10.42,  -99.43, 141.40,  -98.01},
			{114.32, -121.47, 118.14,   82.88, -150.05,  -83.81,  23.35,  -85.82},
			{117.16,  -95.41, 140.40,  -59.35,  -29.23,  -72.39, -25.08,  -76.16},
			{139.20,  -55.96, -32.70,  -68.51,  -26.09,  -74.44, -22.60,  -71.74},
			{-39.62,  -64.73, -39.52,  -65.54,  -38.88,  -66.89, -37.76,  -70.19},
			{-35.34,  -65.03, -38.12,  -66.34,  -29.51,  -89.10,  -2.91,   77.90},
			{-45.29,  -67.44, -27.72,  -87.27,    5.13,   77.49,  30.71,  -93.23},
			{-27.09,  -86.14,   0.30,   59.85,   21.51,  -96.30, 132.67,  -92.91}
		};
	String blockName = "abcdefghijklmnop";
	
	public StructuralAlphabetFingerprint() {
	}
	
	/**
     * Returns a fingerprint as sequence for the given chain. 
     * @param coords
     * @return
     */
	@Override
	public StructuralAlphabetFeature getFingerprint(Point3d[] coordinates) {
		List<Double> angles = new ArrayList<Double>();
		// compute psi and phi angles
		for (int i = 0; i < coordinates.length/3-1; i++) {
			Atom aCA = point3dToAtom(coordinates[i*3]);
			Atom aN = point3dToAtom(coordinates[i*3+1]);
			Atom aC = point3dToAtom(coordinates[i*3+2]);
			Atom bCA = point3dToAtom(coordinates[i*3+3]);
			Atom bN = point3dToAtom(coordinates[i*3+4]);
			Atom bC = point3dToAtom(coordinates[i*3+5]);
			if (aCA == null || aN == null || aC == null || bCA == null || bN == null || bC == null) {
				continue;
			}
			double psi = Calc.torsionAngle(aN,aCA,aC,bN);
			double phi = Calc.torsionAngle(aC,bN,bCA,bC);
			angles.add(psi);
			angles.add(phi);
		}
		String[] assign = new String[(angles.size()-5)/2];
		// compare with the reference blocks
		for (int i = 0; i < angles.size()-7; i+= 2) {
			double[] block = new double[8];
			for (int j = 0; j < 8; j++) {
				block[j] = angles.get(i+j);
			}
			assign[i/2] = assign(block);
		}
		return new StructuralAlphabetFeature(assign);
	}
	
	/**
	 * transfer Point3d to Atom
	 * @param p
	 * @return
	 */
	public Atom point3dToAtom(Point3d p) {
		Atom a = new AtomImpl();
		a.setX(p.x);
		a.setY(p.y);
		a.setZ(p.z);
		return a;
	}
	
	/**
	 * Assign a reference block name to the input block
	 * @param block
	 * @return
	 */
	public String assign(double[] block) {
		double minRmsda = Double.MAX_VALUE;
		char minBlock = 0;
		for (int i = 0; i < references.length; i++) {
			double rmsda = rmsda(references[i],block);
			if (rmsda < minRmsda) {
				minRmsda = rmsda;
				minBlock = blockName.charAt(i);
			}
		}
		return Character.toString(minBlock);
	}

	/**
	 * check rmsd for the angle
	 * @param ref
	 * @param block
	 * @return
	 */
	public double rmsda(double[] ref, double[] block) {
		double sum = 0;
		for (int i = 0; i < ref.length; i++) {
			double dif = angleModule360(ref[i] - block[i]);
			sum += dif * dif;
		}
		return sum;
	}
	
	/**
	 * make sure the angle is between -180 and 180
	 * @param angle
	 * @return
	 */
	public double angleModule360(double angle) {
		if (angle > 180)
			return angle - 360;
		else if (angle < -180)
			return angle + 360;
		else
			return angle;
	}
	
	@Override
	public String getName() {
		return this.getClass().getSimpleName();
	}
}
