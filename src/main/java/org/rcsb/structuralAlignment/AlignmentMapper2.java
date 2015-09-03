package org.rcsb.structuralAlignment;

import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data
 * to a vector of alignment scores
 * 
 * @author  Peter Rose
 */
public class AlignmentMapper2 implements PairFunction<Tuple2<Integer,Integer>,String,Float> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, Point3d[]>>> data = null;

	public AlignmentMapper2(Broadcast<List<Tuple2<String,Point3d[]>>> data) {
		this.data = data;
	}

	/**
	 * Returns a chainId pair with the TM scores
	 */
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Point3d[]> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Point3d[]> t2 = this.data.getValue().get(tuple._2);

		Point3d[] points1 = t1._2;
		Point3d[] points2 = t2._2;
		
		double maxfTm = 0;
		
		int minI = 0;
		int minJ = 0;
		int len = 0;
		
		int inc1 = 1;
		int inc2 = 1;
		int minLength = 32;

		int minLen = minLength/inc2;
		
		// check for null values
		for (int i = 0; i < points1.length;i++) {
			if (points1[i] == null) {
				System.out.println("Null pointer: " + t1._1);
				return new Tuple2<String,Float>(t1._1, 0.0f);
			}
		}
		for (int i = 0; i < points2.length;i++) {
			if (points2[i] == null) {
				System.out.println("Null pointer: " + t2._1);
				return new Tuple2<String,Float>(t2._1, 0.0f);
			}
		}

		for (int i = 0; i < points1.length-minLength;i+=inc1) {
			for (int j = 0; j < points2.length-minLength; j+=inc1) {
				int m = ((Math.min(points1.length-i, points2.length-j)-1)/inc2);

				for (int k = minLen; k < m; k++) {
					Point3d[] p1 = new Point3d[k];
					Point3d[] p2 = new Point3d[k];
					for (int l = 0; l < k; l++) {
						p1[l] = points1[i+l*inc2];
						p2[l] = points2[j+l*inc2];
					}

					SuperPositionQCP qcp = new SuperPositionQCP();
					qcp.set(p1, p2);
					double rmsd = qcp.getRmsd();

					if (rmsd > 5) {
						break;
					}

					double fTm = getfTmScore(k*inc2, rmsd); // calculate approximate TM score
					if (fTm > maxfTm) {
						minI = i;
						minJ = j;
						len = k*inc2;
						maxfTm = fTm;	
					}
				}
			}
		}
		
		double rmsd = 999;
		double tm = 0;
		
		if (maxfTm > 0) {
			
			// calculate exact TM score
			Point3d[] start1 = Arrays.copyOfRange(points1, minI, minI+len);
			Point3d[] start2 = Arrays.copyOfRange(points2, minJ, minJ+len);
			SuperPositionQCP qcp = new SuperPositionQCP();
			qcp.set(start1,  start2);
			rmsd = qcp.getRmsd();
			Point3d[] y = qcp.getTransformedCoordinates();
			tm = SuperPositionQCP.TMScore(start1, y);
		}
	
		String key = t1._1 + minI + "_" +t2._1 + minJ+"_" + len + "_" +rmsd;
//		System.out.println(key + ": " + minI + "-" +maxI + ":" + minJ +"-" + maxJ + " rmsd: " + minRmsd);

		return new Tuple2<String,Float>(key, new Float(tm));
	}
	
	private static double getfTmScore(int alignLen, double rmsd) {
		double d0 = 1.24 * Math.cbrt(alignLen - 15) - 1.8;
		double fTm = alignLen / (alignLen * (1 + (rmsd/d0)*(rmsd/d0)));
		return fTm;
		
	}


}
