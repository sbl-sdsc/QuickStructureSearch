package org.rcsb.structuralAlignment;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.DenseVector;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.apache.spark.mllib.regression.LabeledPoint;
import org.rcsb.fingerprints.GenericFingerprint;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data
 * to a vector of alignment scores
 * 
 * @author  Peter Rose
 */
public class ChainPairToRmsdMapper1 implements Function<Tuple2<Integer,Integer>,LabeledPoint> {
	private static final long serialVersionUID = 1L;
	private Random r = new Random(1);
	private SuperPositionQCP superposer = new SuperPositionQCP();
	private Broadcast<List<Tuple2<String, Point3d[]>>> data = null;
	private GenericFingerprint fingerprint;
	private int fragmentLength;

	public ChainPairToRmsdMapper1(Broadcast<List<Tuple2<String,Point3d[]>>> data, GenericFingerprint fingerprint, int fragmentLength) {
		this.data = data;
		this.fingerprint =fingerprint;
		this.fragmentLength = fragmentLength ;
	}

	/**
	 * Returns a chainId pair with the TM scores
	 */
	public LabeledPoint call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Point3d[]> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Point3d[]> t2 = this.data.getValue().get(tuple._2);
		
		Point3d[] points1 = t1._2;
		Point3d[] points2 = t2._2;
		
		Point3d[] fragment1 = getRandomFragment(points1);
		Point3d[] fragment2 = getRandomFragment(points2);

		double[] hashCodes1 = fingerprint.getFingerprint(fragment1);
		double[] hashCodes2 = fingerprint.getFingerprint(fragment2);
		int index1 = getHashIndex(hashCodes1);
		int index2 = getHashIndex(hashCodes2);

		int match = index1 == index2 ? 1 : 0;
        
		superposer.set(fragment1, fragment2);
		long start = System.nanoTime();
        double feature[] = new double[1];
 //       feature[0] = SuperPosition.calcRmsd(fragment1, fragment2);
        feature[0] = superposer.getRmsd();

        return new LabeledPoint(match,Vectors.dense(feature));
    }
	
	private int getHashIndex(double[] hashCodes) {
		for (int i = 0; i < hashCodes.length; i++) {
			if (hashCodes[i] > 0) {
				return i;
			}
		}
		return -1;
	}
	
	private Point3d[] getRandomFragment(Point3d[] points) {
		int index = getRandomFragmentIndex(points);
		return Arrays.copyOfRange(points, index, index+this.fragmentLength);
	}
	      	
	private int getRandomFragmentIndex(Point3d[] points) {
		int index = r.nextInt(points.length-this.fragmentLength);
		while (hasGaps(points, index)) {
		   index =  r.nextInt(points.length-this.fragmentLength);
		}
		return index;
	}
	/**
	 * Returns true if there is a gap between the C alpha atoms
	 * within a fragment. Note, the first position is not checked, 
	 * since the gap info is for a gap between residue i and i + 1.
	 * @param gaps true if there is a gap in the fragment
	 * @param index start residue of fragment
	 * @return
	 */
	private boolean hasGaps(Point3d[] coords, int index) {
		for (int i = index; i < index+this.fragmentLength; i++) {
			if (coords[i] == null) {
				return true;
			}
		}
		return false;
	}
}
