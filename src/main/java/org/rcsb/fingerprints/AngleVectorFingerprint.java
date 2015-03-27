package org.rcsb.fingerprints;

import java.io.Serializable;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
/**
 * This class generate fingerprint that is an array of angles which is the angle 
 * between each three points
 */
public class AngleVectorFingerprint implements GenericFingerprint, Serializable {

	private static final long serialVersionUID = 1L;

	/**
     * Returns a fingerprint for the given chain. 
     * @param coords coordinates of a macromolecule fragment
     * @return fingerprint
     */
	public double[] getFingerprint(Point3d[] coordinates) {
		// the angle feature array
		double[] features = new double[coordinates.length-1];
		// each two points generate a vector
		Vector3d v1 = new Vector3d();
		Vector3d v2 = new Vector3d();
		for (int i = 1; i < coordinates.length-1; i++) {
			if (coordinates[i-1] == null || coordinates[i] == null || coordinates[i+1] == null)
				continue;
			v1.set(coordinates[i-1]);
			v1.sub(coordinates[i]);
			v2.set(coordinates[i+1]);
			v2.sub(coordinates[i]);
			// the angle betwean two vectors
			features[i] = v2.angle(v1);
		}
		return features;
	}

	public String getName() {
		return this.getClass().getSimpleName();
	}
}
