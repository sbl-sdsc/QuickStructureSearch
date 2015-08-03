package org.rcsb.project3;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the longest local common subsequence over the length of the chains
 * using SmithWaterman algorithm and Gotoh's improvement
 * 
 * @author Chris Li
 */
public class SmithWatermanWithGeoComp implements AlignmentAlgorithmInterface {

	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data = null;
	private Broadcast<List<Tuple2<String,Point3d[]>>> coords = null;
    // print traceback if it is greater than 0
    /* With different open and extend penalty, this class could function the same as LCS or SmithWaterman
     * LCS: open = extend = 0;
     * SmithWaterman = open = extend = 1;
     */
    // open gap penalty
    private double open = 5;
    // extend gap penalty
    private double extend = 0.5;
    private int rotateTime = 4;

    public SmithWatermanWithGeoComp() {
	}
    
    public SmithWatermanWithGeoComp(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
    
    /***
     * Constructor with setting options
     * @param data
     * @param open
     * @param extend
     */
    public SmithWatermanWithGeoComp(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data, double open, double extend) {
		this.data = data;
		this.open = open;
		this.extend = extend;
	}
    
	@Override
	public void setSequence(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
	
	@Override
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) {
		Tuple2<String,SequenceFeatureInterface<?>> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,SequenceFeatureInterface<?>> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		SequenceFeatureInterface<?> v1 = t1._2;
		SequenceFeatureInterface<?> v2 = t2._2;
		
		Alignment<?> SWAlignment = getAlignment(v1, v2, open, extend);
		
		// TODO for OneAgainstAll
		Point3d[] c1 = this.coords.getValue().get(tuple._1)._2;
		Point3d[] c2p = this.coords.getValue().get(tuple._2)._2;
		Point3d[] c2 = new Point3d[c2p.length];
		for (int i = 0; i < c2.length; i++) {
			c2[i] = new Point3d(c2p[i]);
		}
//		Point3d[] c1 = v1.getCoords();
//		Point3d[] c2 = v2.getCoords();
		
		int Lmin = Math.min(c1.length, c2.length);
		
		float value = (float) SWAlignment.calculateScore()/Lmin;

		for (int trial = 0; trial < rotateTime; trial++) {
			SWAlignment = rotateAlign(SWAlignment, c1, c2);
		}

		if (SWAlignment != null) {
			Integer[] v1Order = SWAlignment.getSequence1();
			Integer[] v2Order = SWAlignment.getSequence2();
			
			List<Point3d> lp1 = new ArrayList<Point3d>();
			List<Point3d> lp2 = new ArrayList<Point3d>();
			
			for (int i = 0; i < v1Order.length; i++) {
				if (v1Order[i] != null && v2Order[i] != null) {
					lp1.add(c1[v1Order[i]]);
					lp2.add(c2[v2Order[i]]);
				}
			}
	
			int Laln = lp1.size();
	
			double d0 = 1.24 * Math.cbrt(Lmin - 15.) - 1.8;
			double d0sq = d0 * d0;
	
			double sum = 0;
			for (int i = 0; i < Laln; i++) {
				double d = getDistance(lp1.get(i), lp2.get(i));
				sum += 1./(1 + d * d / d0sq);
			}
			value = (float) (sum/Lmin);
		}
		
		return new Tuple2<String, Float>(key.toString(), (float) value);
	}
	
	/**
	 * Get alignment for the two sequence. Object class casting.
	 * @param v1
	 * @param v2
	 * @param o
	 * @param e
	 * @return
	 */
	@SuppressWarnings("unchecked")
	private <T,K> Alignment<T> getAlignment(SequenceFeatureInterface<T> v1,SequenceFeatureInterface<K> v2,double o, double e) {
		return SmithWatermanGotoh.align(v1, (SequenceFeatureInterface<T>)v2, o, e);
	}
	
	public static double getDistance(Point3d a, Point3d b) {
		double x = a.x - b.x;
		double y = a.y - b.y;
		double z = a.z - b.z;

		double s  = x * x  + y * y + z * z;

		return Math.sqrt(s);
	}
	
	private Alignment<?> rotateAlign(Alignment<?> SWAlignment, Point3d[] c1, Point3d[] c2) {
		if (SWAlignment != null && SWAlignment.getSequence1().length > 0) {
			Integer[] v1Order = SWAlignment.getSequence1();
			Integer[] v2Order = SWAlignment.getSequence2();
			ArrayList<Point3d> lp1 = new ArrayList<Point3d>();
			ArrayList<Point3d> lp2 = new ArrayList<Point3d>();
			for (int i = 0; i < v1Order.length; i++) {
				if (v1Order[i] != null && v2Order[i] != null) {
					lp1.add(c1[v1Order[i]]);
					lp2.add(c2[v2Order[i]]);
				}
			}
			
			Point3d[] p1 = new Point3d[lp1.size()];
			Point3d[] p2 = new Point3d[lp2.size()];
			for (int i = 0; i < p1.length; i++) {
				p1[i] = lp1.get(i);
				p2[i] = lp2.get(i);
			}
			
			SuperPositionQCP qcp = new SuperPositionQCP();
			qcp.set(p1, p2);
			Matrix4d m = qcp.getTransformationMatrix();
			SuperPositionQCP.transform(m, c2);
			
			SequenceFeatureInterface<Point3d> s1 = new Point3dFeature(c1);
			SequenceFeatureInterface<Point3d> s2 = new Point3dFeature(c2);
			return getAlignment(s1, s2, open, extend);
		}
		else 
			return null;
	}
	
	@Override
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> coords) {
		this.coords = coords;
	}
}
