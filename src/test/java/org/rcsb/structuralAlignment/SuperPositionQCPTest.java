package org.rcsb.structuralAlignment;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class SuperPositionQCPTest {
    private double expectedRmsd;
	private Point3d[] x;
    private Point3d[] y;
    private Matrix3d rotation;

	
    /**
     * Test data are taken from: http://theobald.brandeis.edu/qcp/main.c
     */
	@Before
	public void setUp() {
		expectedRmsd = 0.719106;
		
	     Point3d[] f1 = {
			new Point3d(-2.803, -15.373, 24.556),
			new Point3d( 0.893, -16.062, 25.147),
			new Point3d( 1.368, -12.371, 25.885),
			new Point3d(-1.651, -12.153, 28.177),
			new Point3d(-0.440, -15.218, 30.068),
			new Point3d( 2.551, -13.273, 31.372),
			new Point3d( 0.105, -11.330, 33.567)
	     };
	     x = f1;
	     
	     Point3d[] f2 = {
			new Point3d(-14.739, -18.673, 15.040),
			new Point3d(-12.473, -15.810, 16.074),
			new Point3d(-14.802, -13.307, 14.408),
			new Point3d(-17.782, -14.852, 16.171),
			new Point3d(-16.124, -14.617, 19.584),
			new Point3d(-15.029, -11.037, 18.902),
			new Point3d(-18.577, -10.001, 17.996)
			}; 
	     y = f2;
	     
	     rotation = new Matrix3d();
	    	rotation.m00 =  0.72216358; rotation.m01 = 0.69118937; rotation.m02 = -0.02714790;
	        rotation.m10 = -0.52038257; rotation.m11 = 0.51700833; rotation.m12 = -0.67963547;
	        rotation.m20 = -0.45572112; rotation.m21 = 0.50493528; rotation.m22 =  0.73304748;
	}

	/**
	 * Tests the RMSD only calculation
	 */
	@Test
	public void testGetRmsd() {
		SuperPositionQCP superposer = new SuperPositionQCP();
    	superposer.set(x, y);
        double rmsd = superposer.getRmsd();
        assertEquals(expectedRmsd, rmsd, 0.000001);
	}
	
	/**
	 * This test checks the rotation component
	 * of the transformation matrix
	 */
	@Test
	public void testGetTransformation() {
		SuperPositionQCP superposer = new SuperPositionQCP();
    	superposer.set(x, y);
        Matrix4d transformation = superposer.getTransformationMatrix();
        Matrix3d rot = new Matrix3d();
        transformation.getRotationScale(rot);
        assertTrue(rotation.epsilonEquals(rot, 0.00000001));
	}
	
	/**
	 * This tests checks the transformed coordinates by
	 * comparing them with the expected RMSD value.
	 * 
	 * NOTE: Do not calculate the RMSD this way! For the fast calculation
	 * of the RMSD use the getRmsd() method (see testGetRmsd).
	 */
	@Test
	public void testGetTransformedCoordinates() {
		SuperPositionQCP superposer = new SuperPositionQCP();
    	superposer.set(x, y);
        Point3d[] yt = superposer.getTransformedCoordinates();
        double rmsd = SuperPositionQCP.rmsd(x, yt);
        assertEquals(expectedRmsd, rmsd, 0.000001);
	}
}
