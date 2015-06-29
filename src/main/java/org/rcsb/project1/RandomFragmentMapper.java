package org.rcsb.project1;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.structuralAlignment.DistanceRmsd;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/**
 * This class generates two random polymer fragments and calculates the cRMSD and dRMSD.
 * In addition, it calculates the end to end length of the two fragments, the length difference, 
 * and the time to calculate the cRMSD and dRMSD in nanoseconds.
 * 
 * @author Peter Rose
 * @author Dane Malangone
 * @author Justin Li
 * @author Reyd Nguyen
 * @author Joe Sun
 */
public class RandomFragmentMapper implements PairFunction<Tuple2<Integer,Integer>,String,Double[]> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, Point3d[]>>> data = null;
	private int length;
	private int seed;
	private SuperPositionQCP qcp = new SuperPositionQCP();


	public RandomFragmentMapper(Broadcast<List<Tuple2<String,Point3d[]>>> data, int length, int seed) {
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
		
		// Find the last possible start position for a random fragment
		int last1 = points1.length - length;
		int last2 = points2.length - length;

		// Get random fragment start positions (0 .. last-1)
		Random r = new Random(seed);
		int start1 = r.nextInt(last1);
		int start2 = r.nextInt(last2);
		
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

		Double[] rmsds = new Double[7];
		
		// create two random fragments
		Point3d[] fragment1 = Arrays.copyOfRange(points1, start1, start1+length);		
		Point3d[] fragment2 = Arrays.copyOfRange(points2, start2, start2+length);

		// Find the cRMSD
		long t1 = System.nanoTime();
		qcp.set(fragment1, fragment2);
		double crmsd = qcp.getRmsd();	
		long t2 = System.nanoTime();
        rmsds[0] = crmsd;
        
        // Find the dRMSD
        long t3 = System.nanoTime();
        double drmsd = DistanceRmsd.getDistanceRmsd(fragment1, fragment2);
        long t4 = System.nanoTime();
        rmsds[1] = drmsd;
        
        // Calculate end to end distance for the two random fragments
        rmsds[2] = fragment1[0].distance(fragment1[length-1]);
        rmsds[3] = fragment2[0].distance(fragment2[length-1]);
        
        // Calculate the difference in length between the two fragments
        rmsds[4] = Math.abs(rmsds[2]-rmsds[3]);
        
        // Save timing information for cRMSD and dRMSD calculation
        double t12 = t2 - t1;
        double t34 = t4 - t3;
        rmsds[5] = t12;
        rmsds[6] = t34;
		
		return new Tuple2<String, Double[]>(key.toString(), rmsds);
    }
}
