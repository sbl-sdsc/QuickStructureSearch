package org.rcsb.project1;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.structuralAlignment.AfpDistanceRmsd;
import org.rcsb.structuralAlignment.DistanceRmsd;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/**
 * This class generates two random polymer fragments and calculates the cRMSD and dRMSD.
 * In addition, it calculates the end to end length of the two fragments, the length difference, 
 * and the time to calculate the cRMSD and dRMSD in nanoseconds.
 * 
 * @author Peter Rose
 */
public class RandomAfpMapper implements PairFunction<Tuple2<Integer,Integer>,String,Double[]> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, Point3d[]>>> data = null;
	private int length;
	private int seed;
	private static SuperPositionQCP qcp = new SuperPositionQCP();

	public RandomAfpMapper(Broadcast<List<Tuple2<String,Point3d[]>>> data, int length, int seed) {
		this.data = data;
		this.length = length;
		this.seed = seed;
	}

	/**
	 * Returns a tuple2 made of a string and a double array.
	 * 
	 * The string returned consists of the PDBid, ChainId and starting points for both fragments.
	 * The double array consists of the cRMSD and dRMSD values with the time taken for each.
	 */
	public Tuple2<String, Double[]> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Point3d[]> tuple1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Point3d[]> tuple2 = this.data.getValue().get(tuple._2);

		Point3d[] points1 = tuple1._2;
		Point3d[] points2 = tuple2._2;
		
		int maxGapLen = 30;

		// Find the last possible start position for a random fragment
		int last1 = points1.length - length - maxGapLen;
		int last2 = points2.length - length - maxGapLen;
		
		last1 = Math.max(1, last1);
		last2 = Math.max(1, last2);
		
		if (nullPointCoordinateChecker(points1)) {
			throw new IllegalArgumentException("points1: null coordinates");
		}
		if (nullPointCoordinateChecker(points2)) {
			throw new IllegalArgumentException("points2: null coordinates");
		}

		// Get random fragment start positions (0 .. last-1)
		Random r = new Random(seed);
		int start1 = r.nextInt(last1);
		int start2 = start1 + Math.min(8, r.nextInt(maxGapLen)); // 
		int start3 = r.nextInt(last2);
		int start4 = start3 + Math.min(8, r.nextInt(maxGapLen));

		// Create a key consisting of the PdbId.ChainId and 
		// fragment start position of the two polymer chains
		StringBuilder key = new StringBuilder();
		key.append(tuple1._1);
		key.append(",");
		key.append(tuple2._1);
		key.append(",");
		key.append(start1);
		key.append(",");
		key.append(start2);

		Double[] data = new Double[5];

		// create two random aligned fragment pairs
		Point3d[] afp1 = new Point3d[length+length];
		System.arraycopy(points1, start1, afp1, 0, length);
		System.arraycopy(points1, start2, afp1, length, length);
		Point3d[] afp2 = new Point3d[length+length];
		System.arraycopy(points2, start3, afp2, 0, length);
		System.arraycopy(points2, start4, afp2, length, length);

		// Find the cRMSD
		long t1 = System.nanoTime();
		qcp.set(afp1, afp2);
		double crmsd = qcp.getRmsd();
		long t2 = System.nanoTime();
		data[0] = crmsd;

		// create two random fragments
		Point3d[] x1 = Arrays.copyOfRange(points1, start1, start1+length);		
		Point3d[] y1 = Arrays.copyOfRange(points1, start2, start2+length);
		Point3d[] x2 = Arrays.copyOfRange(points2, start3, start3+length);		
		Point3d[] y2 = Arrays.copyOfRange(points2, start4, start4+length);

		// Find the dRMSD
		long t3 = System.nanoTime();
		double drmsd = AfpDistanceRmsd.getAfpDistanceRmsd(x1, y1, x2, y2);
		long t4 = System.nanoTime();
		data[1] = drmsd;

		data[2] = Math.abs(drmsd - crmsd);
		
		// Save timing information for RMSD calculation
		double t12 = t2 - t1;
		double t34 = t4 - t3;
		data[3] = t12;
		data[4] = t34;

		return new Tuple2<String, Double[]>(key.toString(), data);
	}
	
	private static boolean nullPointCoordinateChecker(Point3d[] points) {
		for (Point3d p : points) {
			if (p == null) {
				return true;
			}
		}
		return false;
	}
}
