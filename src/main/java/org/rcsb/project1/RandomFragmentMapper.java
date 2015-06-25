package org.rcsb.project1;

import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.structuralAlignment.DistanceRmsd;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/**
 * This class generates two random fragments and calculates the cRMSD and dRMSD for each pair.
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
	private String fileName;


	public RandomFragmentMapper(Broadcast<List<Tuple2<String,Point3d[]>>> data, int length, int seed, String fileName) {
		this.data = data;
		this.length = length;
		this.seed = seed;
		this.fileName = fileName;
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
		
		int len1 = points1.length - length;
		int len2 = points2.length - length;
		Random r = new Random(seed);
		
		int start1 = r.nextInt(len1);
		int start2 = r.nextInt(len2);
		
		// Make the comma seperated values
		StringBuilder key = new StringBuilder();
		key.append(tuple1._1);
		key.append(",");
		key.append(tuple2._1);
		key.append(",");
		key.append(start1);
		key.append(",");
		key.append(start2);

		Double[] rmsds = new Double[4];

		// New Code
		
		Point3d[] fragment1 = new Point3d[length];
		Point3d[] fragment2 = new Point3d[length];
		
		//Make random (length) long chains
		for (int i = 0; i < length; i++) {
			int place = i + start1;
			fragment1[i] = points1[place];
		}
		
		for (int j = 0; j < length; j++) {
			int place = j + start2;
			fragment2[j] = points2[place];
		}

		// Find the cRMSD
		long t1 = System.nanoTime();
		qcp.set(fragment1, fragment2);
		double crmsd = qcp.getRmsd();
		long t2 = System.nanoTime();
		
        rmsds[0] = crmsd;
        
        // Find the dRMSD
        long t3 = System.nanoTime();
        double drmsd = DistanceRmsd.getDistanceRmsd(fragment1, fragment2);
        rmsds[1] = drmsd;
        long t4 = System.nanoTime();
        
        if (drmsd < 0.5 && crmsd > 1.0) {
    		Point3d[] fragment2Transformed = qcp.getTransformedCoordinates();
        	VisualizeFragmentPair.writeFragmentPair(fragment1, fragment2Transformed, fileName);
        	System.out.println("Found cRMSD vs. dRMSD inconsistency: " + crmsd + ", " + drmsd);
        }
        double t12 = t2 - t1;
        double t34 = t4 - t3;
        rmsds[2] = t12;
        rmsds[3] = t34;
		
		return new Tuple2<String, Double[]>(key.toString(), rmsds);
    }
}
