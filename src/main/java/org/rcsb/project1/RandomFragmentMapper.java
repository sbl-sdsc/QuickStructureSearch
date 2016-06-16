package org.rcsb.project1;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.structuralAlignment.AbsSuperpositionQcp;
import org.rcsb.structuralAlignment.DistanceRmsd;
import org.rcsb.structuralAlignment.IntSuperpositionQcp;
import org.rcsb.structuralAlignment.QCPUpdateable;
import org.rcsb.structuralAlignment.SuperPosition;
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
	private SuperPosition sp = new SuperPosition();
	private QCPUpdateable ucp = new QCPUpdateable();
	private AbsSuperpositionQcp acp = new AbsSuperpositionQcp();
	private IntSuperpositionQcp icp = new IntSuperpositionQcp();


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

		Double[] data = new Double[12];
		
		// create two random fragments
		Point3d[] fragment1 = Arrays.copyOfRange(points1, start1, start1+length);		
		Point3d[] fragment2 = Arrays.copyOfRange(points2, start2, start2+length);

		// Find the cRMSD
		long t1 = System.nanoTime();
		qcp.set(fragment1, fragment2);
		double crmsd = qcp.getRmsd();
		long t2 = System.nanoTime();
        data[0] = crmsd;
        
        // Find rmsd using alternative method
        int[] if1 = getIntCoords(fragment1);
        int[] if2 = getIntCoords(fragment2);
        
        long t7 = System.nanoTime();
 //       double srmsd = sp.calcRmsd(fragment1, fragment2);
        icp.set(if1, if2);
		double srmsd = icp.getRmsd();
        long t8 = System.nanoTime();
        data[1] = srmsd;
        
        // Find the dRMSD
        long t3 = System.nanoTime();
        double drmsd = DistanceRmsd.getDistanceRmsd(fragment1, fragment2);
        long t4 = System.nanoTime();
        data[2] = drmsd;
        
        // Calculate approximate adRMSD
        long t5 = System.nanoTime();
      	double adrmsd = DistanceRmsd.getDistanceRmsdApproximation(fragment1, fragment2);
        long t6 = System.nanoTime();
        data[3] = adrmsd;  
        
        data[4] = Math.abs(drmsd - crmsd);
        // Save timing information for RMSD calculation
        double t12 = t2 - t1;
        double t34 = t4 - t3;
        double t56 = t6 - t5;
        double t78 = t8 - t7;
        data[5] = t12;
        data[6] = t34;
        data[7] = t56;
        data[8] = t78;
        
        
        // Calculate end to end distance for the two random fragments
        data[9] = fragment1[0].distance(fragment1[length-1]);
        data[10] = fragment2[0].distance(fragment2[length-1]);
        
        // Calculate the difference in length between the two fragments
        data[11] = Math.abs(data[9]-data[10]);
		
		return new Tuple2<String, Double[]>(key.toString(), data);
    }
	
	private static int[] getIntCoords(Point3d[] x) {
		int[] iCoord = new int[x.length*3];
		for (int i = 0, c = 0;i < x.length; i++) {
			iCoord[c++] = (int) Math.round(x[i].x *1000);
			iCoord[c++] = (int) Math.round(x[i].y *1000);
			iCoord[c++] = (int) Math.round(x[i].z *1000);
		}
		return iCoord;
	}
}
