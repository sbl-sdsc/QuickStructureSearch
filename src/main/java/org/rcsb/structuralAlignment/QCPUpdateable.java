package org.rcsb.structuralAlignment;

import java.io.Serializable;
import java.util.Arrays;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.apache.commons.math.stat.correlation.Covariance;
import org.apache.commons.math3.stat.correlation.StorelessCovariance;


/**
 * Calculates rmsd and 4x4 transformation matrix for two sets of points.
 * The least-squares rotation is calculated using a quaternion-based 
 * characteristic polynomial (QCP) and a cofactor matrix. The QCP method is
 * currently the fasted known method to calculate the rmsd and superposition
 * of two sets of points.
 * 
 * Usage:
 * 
 * The input consists of 2 Point3d arrays of equal length. The input coordinates
 * are not changed.
 * 
 *    Point3d[] x = ...
 *    Point3d[] y = ...
 *    SuperPositionQCP qcp = new SuperPositionQCP();
 *    qcp.set(x, y);
 *    
 * or with weighting factors [0 - 1]]
 *    double[] weights = ...
 *    qcp.set(x, y, weights);
 *    
 * For maximum efficiency, create a SuperPositionQCP object once and reuse it.
 * 
 * A. Calculate rmsd only
 * 	  double rmsd = qcp.getRmsd();
 * 
 * B. Calculate a 4x4 transformation (rotation and translation) matrix
 *    Matrix4d rottrans = qcp.getTransformationMatrix();
 * 
 * C. Get transformated points (y superposed onto the reference x)
 *    Point3d[] ySuperposed = qcp.getTransformedCoordinates();
 * 
 * 
 * Citations
 * 
 * Liu P, Agrafiotis DK, & Theobald DL (2011)
 * Reply to comment on: "Fast determination of the optimal rotation matrix for macromolecular superpositions."
 * Journal of Computational Chemistry 32(1):185-186. [http://dx.doi.org/10.1002/jcc.21606]
 *
 * Liu P, Agrafiotis DK, & Theobald DL (2010)
 * "Fast determination of the optimal rotation matrix for macromolecular superpositions."
 * Journal of Computational Chemistry 31(7):1561-1563. [http://dx.doi.org/10.1002/jcc.21439]
 *
 * Douglas L Theobald (2005)
 * "Rapid calculation of RMSDs using a quaternion-based characteristic polynomial."
 * Acta Crystallogr A 61(4):478-480. [http://dx.doi.org/10.1107/S0108767305015266 ]
 * 
 * This is adoption of the original C code QCProt 1.4 (2012, October 10) to Java. 
 * The original C source code is available from http://theobald.brandeis.edu/qcp/ and was developed by
 * 
 * Douglas L. Theobald
 * Department of Biochemistry
 * MS 009
 * Brandeis University
 * 415 South St
 * Waltham, MA  02453
 * USA
 *
 * dtheobald@brandeis.edu
 *                  
 * Pu Liu
 * Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 * 665 Stockton Drive
 * Exton, PA  19341
 * USA
 *
 * pliu24@its.jnj.com
 * 
 * The following notice is from the original C source code:
 * 
 * Copyright (c) 2009-2013 Pu Liu and Douglas L. Theobald
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are permitted
 * provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list
 *    of conditions and the following disclaimer in the documentation and/or other materials
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
 *    endorse or promote products derived from this software without specific prior written
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 * @author Peter Rose (adopted to Java)
 */
public final class QCPUpdateable implements Serializable {
	private static final long serialVersionUID = 1L;
	private static final double EVE_PREC = 1E-6;
    private static final double EVAL_PREC = 1E-11;
    
    private Point3d[] x;
    private Point3d[] y;
    private double[] weight;
    
    private Point3d xCentroid;
    private Point3d yCentroid; 
    private Matrix4d transformation = new Matrix4d();
    private double rmsd = 0;
    private double upperBound;
    private double g;
    private double Sxy, Sxz, Syx, Syz, Szx, Szy;
    private double SxxpSyy, Szz, mxEigenV, SyzmSzy,SxzmSzx, SxymSyx;
    private double SxxmSyy, SxypSyx, SxzpSzx;
    private double Syy, Sxx, SyzpSzy;
    private double Cxx, Cxy, Cxz, Cyx, Cyy, Cyz, Czx, Czy, Czz, Cg; 
    private boolean rmsdCalculated = false;
    private boolean transformationCalculated = false;
    private boolean centered = false;
    private int length;
    
    /**
     * Default constructor
     */
    public QCPUpdateable() {
    	this.centered = false;
    }
    /**
     * Constructor with option to set centered flag. This constructor
     * should be used if both coordinate input set have been centered at the origin.
     * @param centered if set true, the input coordinates are already centered at the origin
     */
    public QCPUpdateable(boolean centered) {
		this.centered = centered;
	}
	/**
     * Sets the two input coordinate arrays. These input arrays must be of
     * equal length. Input coordinates are not modified.
     * @param x 3d points of reference coordinate set
     * @param y 3d points of coordinate set for superposition
     */
    public void set(Point3d[] x, Point3d[] y) {
    	this.x = x;
    	this.y = y;
    	this.weight = null;
    	this.length = x.length;
        rmsdCalculated = false;
        transformationCalculated = false;
    }
    
    /**
     * Sets the two input coordinate arrays and weight array. 
     * All input arrays must be of equal length. 
     * Input coordinates are not modified.
     * @param x 3d points of reference coordinate set
     * @param y 3d points of coordinate set for superposition
     * @param weight a weight in the inclusive range [0,1] for each point
     */
    public void set(Point3d[] x, Point3d[] y, double[] weight) {
    	this.x = x;
    	this.y = y;
    	this.weight = weight;
    	this.length = x.length;
        rmsdCalculated = false;
        transformationCalculated = false;
    }
    
    /**
     * Return the RMSD of the superposition of input coordinate set y onto x.
     * Note, this is the fasted way to calculate an RMSD without actually
     * superposing the two sets. The calculation is performed "lazy", meaning
     * calculations are only performed if necessary.
     * @return root mean square deviation for superposition of y onto x
     */
    public double getRmsd() {
    	if (! rmsdCalculated) {
    		calcRmsd(x, y);
    	}
    	return rmsd;
    }
    
    /**
     * Fast approximate implementation of TM score.
     * Bioinf. (2012) 28, 1209-1215
     * @param rmsd
     * @return
     */
	public double getfTmScore(double rmsd) {
		if (! rmsdCalculated) {
    		calcRmsd(x, y);
    	}
		double d0 = 1.24 * Math.cbrt(length - 15) - 1.8;
		return length / (length * (1 + (rmsd/d0)*(rmsd/d0)));
	}
    
    /**
     * Returns a 4x4 transformation matrix that transforms the y coordinates onto the x coordinates.
     * The calculation is performed "lazy", meaning calculations are only performed if necessary.
     * @return 4x4 transformation matrix to transform y coordinates onto x
     */
    public Matrix4d getTransformationMatrix() {
    	if (! transformationCalculated) {
    		if (! rmsdCalculated) {
    			calcRmsd(x, y);
    		}
    		Matrix3d rotmat = calcRotationMatrix();
    		if (! centered) {
    			calcTransformation(rotmat);
    		} else {
    			transformation.setIdentity();
    			transformation.set(rotmat);
    		}
    		transformationCalculated = true;
    	}
    	return transformation; 	
    }
    
    /**
     * Returns the transformed (superposed) y coordinates
     * @return transformed y coordinates
     */
    public Point3d[] getTransformedCoordinates() {
    	Matrix4d matrix = getTransformationMatrix();
    	Point3d[] points = new Point3d[y.length];
    	for (int i = 0; i < y.length; i++) {
    		points[i] = new Point3d(y[i]);
    		matrix.transform(points[i]);
    	}
    	return points;
    }
    
    /**
     * Calculates the RMSD value for superposition of y onto x.
     * @param x 3d points of reference coordinate set
     * @param y 3d points of coordinate set for superposition
     */
    private void calcRmsd(Point3d[] x, Point3d[] y) {
    	if (centered) {
        	xCentroid = new Point3d(0,0,0);
        	yCentroid = new Point3d(0,0,0);
    	} else {
    		xCentroid = centroid(x);
    		yCentroid = centroid(y);
    	}
    	if (weight == null) {
//    		innerProduct();
    		innerProductAlt();
    	} else {
    		innerProductWeighted();
    	}
    	calcRmsd(x.length);
    	rmsdCalculated = true;
    }
    
    
    /**
     * Calculates the RMSD value for superposition of y onto x.
     * @param x 3d points of reference coordinate set
     * @param y 3d points of coordinate set for superposition
     */
    public void updateRmsd(Point3d px, Point3d py) {
//    	long t1 = System.nanoTime();
    	Sxx -=  Cxx;
    	Sxy -=  Cxy;
    	Sxz -=  Cxz;

    	Syx -=  Cyx;
    	Syy -=  Cyy;
    	Syz -=  Cyz;

    	Szx -=  Czx;
    	Szy -=  Czy;
    	Szz -=  Czz;
    	
    	g -= Cg;
    	
       length++;
       xCentroid = updateCentroid(xCentroid, px);
       yCentroid = updateCentroid(yCentroid, py);
       
    	if (weight == null) {
    		updateInnerProduct(px, py);
    	} 
//    	long t2 = System.nanoTime();
//    	long time = t2 - t1;
  //  	System.out.println("Update inner product: " + time);
    	
 //   	t1 = System.nanoTime();
    	calcRmsd(length);
//    	t2 = System.nanoTime();
//    	time = t2 - t1;
//      	System.out.println("Update rmsd: " + time);
    	rmsdCalculated = true;
    }
    
    private Point3d updateCentroid(Point3d centroid, Point3d p) {
    	Point3d p1 = new Point3d(p);
    	p1.sub(centroid);
        p1.scale(1.0/(length));
    	p1.add(centroid);
    	return p1;
    }
 
    /**
     * Calculates a 4x4 transformation matrix to superpose coordinate set y onto x
     * @param rotmat rotation matrix for superposition
     */
    private void calcTransformation(Matrix3d rotmat) {
    	// set rotation
    	transformation.setIdentity();
        transformation.set(rotmat);

        // combine with y -> origin translation
        Matrix4d trans = new Matrix4d();
        trans.setIdentity();
        Vector3d yv = new Vector3d(yCentroid);
        yv.negate();
        trans.setTranslation(yv);
        transformation.mul(transformation, trans);
        
        // combine with origin -> x translation  
        Matrix4d transInverse = new Matrix4d(); 
        transInverse.setIdentity();    
        Vector3d xv = new Vector3d(xCentroid);
        transInverse.setTranslation(xv);
        transformation.mul(transInverse, transformation);
    }
    
    /** 
     * Calculates the inner product between two coordinate sets x and y. It also
     * calculates an upper bound of the most positive root of the key matrix.
     * @return
     */
    private void innerProduct() {
        g = 0.0;
        
        Sxx = 0;
        Sxy = 0;
        Sxz = 0;
        Syx = 0;
        Syy = 0;
        Syz = 0;
        Szx = 0;
        Szy = 0;
        Szz = 0; 

        for (int i = 0, n = x.length; i < n; i++)
        { 
        	double x1 = x[i].x - xCentroid.x;
        	double y1 = x[i].y - xCentroid.y;
        	double z1 = x[i].z - xCentroid.z;

        	double x2 = y[i].x - yCentroid.x;
        	double y2 = y[i].y - yCentroid.y;
        	double z2 = y[i].z - yCentroid.z;

        	g += x1*x1 + y1*y1 + z1*z1 + x2*x2 + y2*y2 + z2*z2;

        	Sxx +=  x1*x2;
        	Sxy +=  x1*y2;
        	Sxz +=  x1*z2;

        	Syx +=  y1*x2;
        	Syy +=  y1*y2;
        	Syz +=  y1*z2;

        	Szx +=  z1*x2;
        	Szy +=  z1*y2;
        	Szz +=  z1*z2;  
        }

        upperBound = g * 0.5;
    }
    
    /** 
     * Adds two points to the inner product matrix. It also
     * calculates an upper bound of the most positive root of the key matrix.
     * @return
     */
    private void updateInnerProduct(Point3d px, Point3d py) {
    	initializeInterproductMatrix();
    	
    	double x1 = px.x;
    	double y1 = px.y;
    	double z1 = px.z;

    	double x2 = py.x;
    	double y2 = py.y;
    	double z2 = py.z;
    	
        g += Cg + x1*x1 + y1*y1 + z1*z1 + x2*x2 + y2*y2 + z2*z2;

        Sxx += Cxx + x1*x2;
        Sxy += Cxy + x1*y2;
        Sxz += Cxz + x1*z2;

        Syx += Cyx + y1*x2;
        Syy += Cyy + y1*y2;
        Syz += Cyz + y1*z2;

        Szx += Czx + z1*x2;
        Szy += Czy + z1*y2;
        Szz += Czz + z1*z2;  
        
        upperBound = g * 0.5;
    }

    
    /** 
     * Calculates the inner product between two coordinate sets x and y. It also
     * calculates an upper bound of the most positive root of the key matrix.
     * @return
     */
    private void innerProductAlt() {
        initializeInterproductMatrix();
    	Sxx =  Cxx;
    	Sxy =  Cxy;
    	Sxz =  Cxz;

    	Syx =  Cyx;
    	Syy =  Cyy;
    	Syz =  Cyz;

    	Szx =  Czx;
    	Szy =  Czy;
    	Szz =  Czz;
    	
    	g = Cg;
        
        for (int i = 0, n = x.length; i < n; i++)
        { 
        	double x1 = x[i].x;
        	double y1 = x[i].y;
        	double z1 = x[i].z;

        	double x2 = y[i].x;
        	double y2 = y[i].y;
        	double z2 = y[i].z;
 
            g += x1*x1 + y1*y1 + z1*z1 + x2*x2 + y2*y2 + z2*z2;

        	Sxx +=  x1*x2;
        	Sxy +=  x1*y2;
        	Sxz +=  x1*z2;

        	Syx +=  y1*x2;
        	Syy +=  y1*y2;
        	Syz +=  y1*z2;

        	Szx +=  z1*x2;
        	Szy +=  z1*y2;
        	Szz +=  z1*z2;  	
        }

        upperBound = g * 0.5;
    }
    
	private void initializeInterproductMatrix() {
		Cxx = -xCentroid.x*yCentroid.x*length;
    	Cxy = -xCentroid.x*yCentroid.y*length;
    	Cxz = -xCentroid.x*yCentroid.z*length;

    	Cyx = -xCentroid.y*yCentroid.x*length;
    	Cyy = -xCentroid.y*yCentroid.y*length;
    	Cyz = -xCentroid.y*yCentroid.z*length;

    	Czx = -xCentroid.z*yCentroid.x*length;
    	Czy = -xCentroid.z*yCentroid.y*length;
    	Czz = -xCentroid.z*yCentroid.z*length;
    	
    	double xx1 = -xCentroid.x*xCentroid.x;
       	double yy1 = -xCentroid.y*xCentroid.y;
       	double zz1 = -xCentroid.z*xCentroid.z;
       	
    	double xx2 = -yCentroid.x*yCentroid.x;
       	double yy2 = -yCentroid.y*yCentroid.y;
       	double zz2 = -yCentroid.z*yCentroid.z;

        Cg = (xx1 + yy1 + zz1 + xx2 + yy2 + zz2) * length;
	}
    
    /** 
     * Calculates the inner product between two weighted coordinate sets x and y. It also
     * calculates an upper bound of the most positive root of the key matrix.
     */
	private void innerProductWeighted() {
		double g1 = 0.0;
		double g2 = 0.0;

		Sxx = 0;
		Sxy = 0;
		Sxz = 0;
		Syx = 0;
		Syy = 0;
		Syz = 0;
		Szx = 0;
		Szy = 0;
		Szz = 0;

		for (int i = 0, n = x.length; i < n; i++) {
			double x1 = weight[i] * (x[i].x - xCentroid.x);
			double y1 = weight[i] * (x[i].y - xCentroid.y);
			double z1 = weight[i] * (x[i].z - xCentroid.z);

			g1 += x1 * (x[i].x - xCentroid.x) + y1 * (x[i].y - xCentroid.y) + z1 * (x[i].z - xCentroid.z);

			double x2 = y[i].x - yCentroid.x;
			double y2 = y[i].y - yCentroid.y;
			double z2 = y[i].z - yCentroid.z;

			g2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

			Sxx += x1 * x2;
			Sxy += x1 * y2;
			Sxz += x1 * z2;

			Syx += y1 * x2;
			Syy += y1 * y2;
			Syz += y1 * z2;

			Szx += z1 * x2;
			Szy += z1 * y2;
			Szz += z1 * z2;
		}

		upperBound = (g1 + g2) * 0.5;
	}

	public double getUpperBoundRmsd() {
		return Math.sqrt(Math.abs(2.0 * upperBound/length));
	}
	
	/**
	 * Calculates the RMSD value by determining the most positive root of the key matrix
	 * using the Newton-Raphson method.
	 * @param len length of the input coordinate arrays
	 */
    private void calcRmsd(int len) {     
        double Sxx2 = Sxx * Sxx;
        double Syy2 = Syy * Syy;
        double Szz2 = Szz * Szz;

        double Sxy2 = Sxy * Sxy;
        double Syz2 = Syz * Syz;
        double Sxz2 = Sxz * Sxz;

        double Syx2 = Syx * Syx;
        double Szy2 = Szy * Szy;
        double Szx2 = Szx * Szx;

        double SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
        double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

        double c2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
        double c1 = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

        SxzpSzx = Sxz + Szx;
        SyzpSzy = Syz + Szy;
        SxypSyx = Sxy + Syx;
        SyzmSzy = Syz - Szy;
        SxzmSzx = Sxz - Szx;
        SxymSyx = Sxy - Syx;
        SxxpSyy = Sxx + Syy;
        SxxmSyy = Sxx - Syy;
        
        double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

        double c0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

        mxEigenV = upperBound;      
 
        // calculate most positive root with the Newton-Raphson method
        int iter;
        for (iter = 0; iter < 50; iter++)
        {
            double oldg = mxEigenV;
            double x2 = mxEigenV*mxEigenV;
            double b = (x2 + c2)*mxEigenV;
            double a = b + c1;
            double delta = ((a*mxEigenV + c0)/(2.0*x2*mxEigenV + b + a));
            mxEigenV -= delta;
        
            if (Math.abs(mxEigenV - oldg) < Math.abs(EVAL_PREC*mxEigenV)) 
                break;
        }

        if (iter == 50) 
           System.err.println("SuperPositionQCP: Newton-Raphson not converged after " + iter + " iterations");

        // use absolute value to guard against extremely small, 
        // but *negative* numbers due to floating point error
        rmsd = Math.sqrt(Math.abs(2.0 * (upperBound - mxEigenV)/len));
    }
    
    /**
     * Calculates the rotation matrix to superpose y onto x.
     * @return 3x3 rotation matrix
     */
    private Matrix3d calcRotationMatrix() {
        double a11 = SxxpSyy + Szz-mxEigenV;
        double a12 = SyzmSzy; 
        double a13 = - SxzmSzx; 
        double a14 = SxymSyx;
        double a21 = SyzmSzy; 
        double a22 = SxxmSyy - Szz-mxEigenV; 
        double a23 = SxypSyx; 
        double a24 = SxzpSzx;
        double a31 = a13; 
        double a32 = a23; 
        double a33 = Syy-Sxx-Szz - mxEigenV; 
        double a34 = SyzpSzy;
        double a41 = a14; 
        double a42 = a24; 
        double a43 = a34; 
        double a44 = Szz - SxxpSyy - mxEigenV;
        double a3344_4334 = a33 * a44 - a43 * a34; 
        double a3244_4234 = a32 * a44-a42*a34;
        double a3243_4233 = a32 * a43 - a42 * a33; 
        double a3143_4133 = a31 * a43-a41*a33;
        double a3144_4134 = a31 * a44 - a41 * a34; 
        double a3142_4132 = a31 * a42-a41*a32;
        double q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
        double q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
        double q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
        double q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

        double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

       // The following code tries to calculate another column in the adjoint matrix when the norm of the 
       // current column is too small. Usually this block will never be activated.
        if (qsqr < EVE_PREC) {
            q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
            q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
            q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
            q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

            if (qsqr < EVE_PREC) {
                double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
                double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
                double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

                q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
                q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
                q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
                q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

                if (qsqr < EVE_PREC) {
                    q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                    q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                    q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                    q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                    qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                    
                    if (qsqr < EVE_PREC) {
                        // qsqr is still too small, return the identity matrix.
                    	Matrix3d rotmat = new Matrix3d();
                        rotmat.setIdentity();
                    	return rotmat;
                    }
                }
            }
        }
        
        return toRotationMatrix(q1, q2, q3, q4, qsqr);
    }

	private Matrix3d toRotationMatrix(double q1, double q2, double q3, double q4, double qsqr) {
		double normq = Math.sqrt(qsqr);
        q1 /= normq;
        q2 /= normq;
        q3 /= normq;
        q4 /= normq;

        double a2 = q1 * q1;
        double x2 = q2 * q2;
        double y2 = q3 * q3;
        double z2 = q4 * q4;

        double xy = q2 * q3;
        double az = q1 * q4;
        double zx = q4 * q2;
        double ay = q1 * q3;
        double yz = q3 * q4;
        double ax = q1 * q2;

        Matrix3d rotmat = new Matrix3d();
        rotmat.m00 = a2 + x2 - y2 - z2;
        rotmat.m01 = 2 * (xy + az);
        rotmat.m02 = 2 * (zx - ay);

        rotmat.m10 = 2 * (xy - az);
        rotmat.m11 = a2 - x2 + y2 - z2;
        rotmat.m12 = 2 * (yz + ax);

        rotmat.m20 = 2 * (zx + ay);
        rotmat.m21 = 2 * (yz - ax);
        rotmat.m22 = a2 - x2 - y2 + z2;
        
        return rotmat;
	}
	
    public static Point3d centroid(Point3d[] x) {
        Point3d center = new Point3d();
        for (Point3d p: x) {
            center.add(p);
        }
        center.scale(1.0/x.length);
        return center;
    }
    
    public static void center(Point3d[] x) {
        Point3d center = centroid(x);
        center.negate();
        translate(center, x);
    }
    
    public static void translate(Point3d trans, Point3d[] x) {
        for (Point3d p: x) {
            p.add(trans);
        }
    }
    
    public static void transform(Matrix4d rotTrans, Point3d[] x) {
        for (Point3d p: x) {
            rotTrans.transform(p);
        }
    }
    
    public static double rmsd(Point3d[] x, Point3d[] y) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i].distanceSquared(y[i]);
        }
        return (double)Math.sqrt(sum/x.length);
    }
}
