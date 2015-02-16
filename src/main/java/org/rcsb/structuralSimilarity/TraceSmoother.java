package org.rcsb.structuralSimilarity;

import java.io.IOException;
import java.io.PrintWriter;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;

public class TraceSmoother {
	private static AtomCache cache = initializeCache();

	public static void main(String[] args) throws IOException {
		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/peter/Data/PDB/");
		String fileName = "/Users/peter/Data/TraceSmooth/";

		System.out.println("TraceSmoother");
//		String pdbId = "4QTS.B"; // high writhing number	
//		String pdbId = "1OHR.A"; // non-perfect smoothing
//		String pdbId = "1STP.A"; // protein with gaps
		String pdbId = "3LB2.A"; // 
//		String pdbId = "1IRD.A"; // 
//		String pdbId = "1EP4.A"; // 
//		String pdbId = "2YNG.A"; // 

		int iterations = 2;
		Structure s = null;
		try {
			s = StructureIO.getStructure(pdbId);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Chain c = s.getChains().get(0);

		Atom[] ca = StructureTools.getAtomCAArray(c);
		
		
		fileName += pdbId;
		
		// save original coordinates
		PrintWriter writer = new PrintWriter(fileName + "0.pdb");
		Point3d[] points = new Point3d[ca.length];
		for (int i = 0; i < points.length; i++) {
			points[i] = new Point3d(ca[i].getCoords());
			writer.println(ca[i].toPDB());
		}
		writer.close();
		
		// save savitzky golay smoothed coordinates (cubic 7 point, 2 iterations)
		writer = new PrintWriter(fileName + "sg5.pdb");
		Point3d[] sPoints = smoothSavitzkyGolayCubic7Point(points, 2);
		int offset = (points.length - sPoints.length)/2;

		for (int i = offset; i < points.length-offset; i++) {
			sPoints[i-offset].get(ca[i].getCoords());
			writer.println(ca[i].toPDB());
		}
		
		writer.close();
	}
	
	public static Point3d[] smoothRogen(Point3d[] points, int iterations) {		
		double[] x = new double[points.length];
		double[] y = new double[points.length];
		double[] z = new double[points.length];
		
		for (int i = 0; i < points.length; i++) {
			x[i] = points[i].x;
			y[i] = points[i].y;
			z[i] = points[i].z;
		}
		
		for (int i = 0; i < iterations; i++) {
			x = smoothRogen(x);
			y = smoothRogen(y);
			z = smoothRogen(z);
		}

		Point3d[] smoothedPoints = new Point3d[x.length];
		for (int i = 0; i < x.length; i++) {
			smoothedPoints[i] = new Point3d(x[i],y[i],z[i]);
		}
		
		return smoothedPoints;
	}
	
	/**
	 * P. Rogen, Evaluating protein structure descriptors and tuning Gauss integral
	 * based descriptors, J. Phys.: Condens. Matter (2005) 17, S1523-S1538
	 * @param x
	 * @return
	 */
	private static double[] smoothRogen(double[] x) {	    	
		double[] y = new double[x.length];
        double a = 2.4;
        double b = 2.1;
        
		y[0] = x[0];
		y[1] = (a*x[0] + b*x[1] + a*x[2])/(a+a+b);
		for (int i = 2; i < x.length-2; i++) {
			y[i] = (x[i-2]+ a*x[i-1] + b*x[i] + a*x[i+1] + x[i+2])/(2+a+a+b);
		}
		y[x.length-2] =  (a*x[x.length-3] + b*x[x.length-2] + a*x[x.length-1])/(a+a+b);
		y[x.length-1] = x[x.length-1];

		return y;
	}
	
	/**
	 * http://en.wikipedia.org/wiki/Savitzky–Golay_filter
	 * Note, the smoothed curve will have fewer points in its current implementation
	 * @param points
	 * @param iterations
	 * @return
	 */
	public static Point3d[] smoothSavitzkyGolayLinear4Point(Point3d[] points, int iterations) {	
		int[] coefficients = {1,1,1,1};
		int norm = 4;
		return smoothSavitzkyGolay(points, coefficients, norm, iterations);
	}
	
	/**
	 * http://en.wikipedia.org/wiki/Savitzky–Golay_filter
	 * Note, the smoothed curve will have fewer points in its current implementation
	 * @param points
	 * @param iterations
	 * @return
	 */
	public static Point3d[] smoothSavitzkyGolayCubic7Point(Point3d[] points, int iterations) {	
		int[] coefficients = {-2,3,6,7,6,3,-2};
		int norm = 21;
		return smoothSavitzkyGolay(points, coefficients, norm, iterations);
	}
	
	public static Point3d[] smoothSavitzkyGolay(Point3d[] points, int[] coefficients, int norm, int iterations) {	
		double[] x = new double[points.length];
		double[] y = new double[points.length];
		double[] z = new double[points.length];
		
		for (int i = 0; i < points.length; i++) {
			x[i] = points[i].x;
			y[i] = points[i].y;
			z[i] = points[i].z;
		}
		
		for (int i = 0; i < iterations; i++) {
			x = smoothSavitzkyGolay(x, coefficients, norm);
			y = smoothSavitzkyGolay(y, coefficients, norm);
			z = smoothSavitzkyGolay(z, coefficients, norm);
		}

		// there will be fewer points after smoothing
		Point3d[] smoothedPoints = new Point3d[x.length];
		for (int i = 0; i < x.length; i++) {
			smoothedPoints[i] = new Point3d(x[i],y[i],z[i]);
		}
		
		return smoothedPoints;
	}
	
    private static double[] smoothSavitzkyGolay(double[] x, int[] coefficients, int norm) {
    	int n = coefficients.length - 1;
    	int m = x.length - n;
    	
    	double[] y = new double[m];
    	double[] p = new double[coefficients.length];
    	
        for (int i = 1; i < coefficients.length; i++) {
        	p[i] = x[i-1];	
        }
        for (int i = 0; i < m; i++) {
        	for (int k = 0; k < n; k++) {
        		p[k]= p[k+1];
        	}
        	p[n] = x[i+n];
        	double sum = 0;
        	for (int k = 0; k < coefficients.length; k++) {
        		sum += coefficients[k] * p[k];
        	}
        	y[i] = sum/norm;
        }
        return y;
    }

	private static AtomCache initializeCache() {
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
		cache.setFileParsingParams(params);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		return cache;
	}
}
