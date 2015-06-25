package org.rcsb.structuralSimilarity;

import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.mllib.linalg.Vector;

import scala.Tuple2;

/**
 * This class maps a pair of chains, specified by two indices into the broadcasted data
 * to a vector of alignment scores
 * 
 * @author  Peter Rose
 */
public class ChainPairToTmMapper implements PairFunction<Tuple2<Integer,Integer>,String,Float[]> {
	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String, Point3d[]>>> data = null;


	public ChainPairToTmMapper(Broadcast<List<Tuple2<String,Point3d[]>>> data) {
		this.data = data;
	}

	/**
	 * Returns a chainId pair with the TM scores
	 */
	public Tuple2<String, Double[]> call(Tuple2<Integer, Integer> tuple) throws Exception {
		Tuple2<String,Point3d[]> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,Point3d[]> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		Point3d[] points1 = t1._2;
		Point3d[] points2 = t2._2;

		int len1 = points1.length - 8;
		int len2 = points2.length - 8;
		Random r = new Random();
		
		int start1 = r.nextInt(len1);
		int start2 = r.nextInt(len2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		key.append(",");
		key.append(start1);
		key.append(",");
		key.append(start2);
		
		Double[] rmsds = new Double[2];
		
		// rmsds[0] =  . . crmsd
		// rmsds[1] =  . . drmsd
		
		System.out.println(key);
		
		return new Tuple2<String, Double[](key.toString(), rmsds)>;
		
    }
}
