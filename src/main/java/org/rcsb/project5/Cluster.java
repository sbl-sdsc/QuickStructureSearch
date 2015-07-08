package org.rcsb.project5;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

public class Cluster {
	private int seqClusterId;
	private int strClusterId;
	private List<Tuple2<String, SimplePolymerChain>> strCluster;
	private Tuple2<String, SimplePolymerChain> repChain;
	private SuperPositionQCP qcp;
	private LongestCommonSubstring lcs;

	public Cluster(int seqClusterId, int strClusterId, List<Tuple2<String, SimplePolymerChain>> strCluster, Tuple2<String, SimplePolymerChain> repChain) {
		this.seqClusterId = seqClusterId;
		this.strClusterId = strClusterId;
		this.strCluster = strCluster;
		this.repChain = repChain;
		qcp = new SuperPositionQCP();
		lcs = new LongestCommonSubstring();
	}

	public int getSeqClusterId() {
		return seqClusterId;
	}

	public int getStrClusterId() {
		return strClusterId;
	}

	public List<Tuple2<String, SimplePolymerChain>> getStrCluster() {
		return strCluster;
	}

	public Tuple2<String, SimplePolymerChain> getRepChain() {
		return repChain;
	}

	public int size() {
		return strCluster.size();
	}

	public void setRepChain(Tuple2<String, SimplePolymerChain> tuple) {
		repChain = tuple;
	}

	public void findRepChain() {
		double[] scores = new double[size()];
		double tempScore;
		int gaps;
		int holes;
		for(int n = 0; n < size(); n ++) {
			List<Tuple2<String, SimplePolymerChain>> a = getStrCluster();
			Tuple2<String, SimplePolymerChain> b = a.get(n);
			SimplePolymerChain c = b._2;
			Point3d[] d = c.getCoordinates();
			double e = d.length;
			tempScore = getStrCluster().get(n)._2.getCoordinates().length;
			Point3d[] coordinates = strCluster.get(n)._2.getCoordinates();
			gaps = 0;
			holes = 0;
			boolean newGap = true;
			for(int index = 0; index < coordinates.length; index ++) {
				if(coordinates[index] == null) {
					if(newGap)
						gaps ++;
					newGap = false;
				} else
					newGap = true;
			}
			for(int pos = 0; pos < coordinates.length; pos ++) {
				if(coordinates[pos] == null)
					holes ++;
			}
			tempScore = tempScore - gaps - 0.1*holes;
			scores[n] = tempScore;
		}
		double highScore = scores[0];
		List<Integer> highScoreIndices = new ArrayList<Integer>();
		highScoreIndices.add(0);
		for(int count = 1; count < scores.length; count ++) {
			if(scores[count] > highScore) {
				highScore = scores[count];
				highScoreIndices = new ArrayList<Integer>();
				highScoreIndices.add(count);
			}
			else if(scores[count] == highScore) {
				highScoreIndices.add(count);
			}
		}
		if(highScoreIndices.size() == 1) {
			setRepChain(getStrCluster().get(highScoreIndices.get(0)));
		}
		else {
			List<Tuple2<String, SimplePolymerChain>> tiebreakerChains = new ArrayList<Tuple2<String, SimplePolymerChain>>();
			for(int x = 0; x < highScoreIndices.size(); x ++) {
				tiebreakerChains.add(getStrCluster().get(highScoreIndices.get(x)));
			}
			tiebreaker(tiebreakerChains);
		}
	}

	public void tiebreaker(List<Tuple2<String, SimplePolymerChain>> tiebreakerList) {
		double[][] array = new double[tiebreakerList.size()][tiebreakerList.size()];
		List<Integer> startEnd = null;
		double cRmsd;
		for (int outer = 0; outer < array.length - 1; outer++) {
			for (int inner = outer + 1; inner < array.length; inner++) {
				startEnd = lcs.longestCommonSubstring(
						tiebreakerList.get(outer)._2.getSequence(), tiebreakerList.get(inner)._2.getSequence());
				cRmsd = getcRmsd(tiebreakerList.get(outer), tiebreakerList.get(inner), startEnd.get(0),
						startEnd.get(1), startEnd.get(2), startEnd.get(3));
				array[outer][inner] = cRmsd;
			}
		}
		double[] simValue = new double[tiebreakerList.size()];
		for (int outer = 0; outer < simValue.length; outer++) {
			for (int inner = 0; inner < simValue.length; inner++) {
				if (inner < outer) {
					simValue[outer] += array[outer][inner];
				} else if (inner > outer) {
					simValue[outer] += array[inner][outer];
				}
			}
		}

		double min = simValue[0];
		int minIndex = 0;
		for (int i = 1; i < simValue.length; i++) {
			if(simValue[i] < min) {
				min = simValue[i];
				minIndex = i;
			}
		}
		setRepChain(tiebreakerList.get(minIndex));
	}

	public String toString() {
		String list = "";
		if (repChain != null) {
			list += seqClusterId + "." + strClusterId + "(" + repChain._1 + "): ";
		}else {
			list += seqClusterId + "." + strClusterId + "(null): ";
		}
		for(int i = 0; i < size(); i++) {
			list += strCluster.get(i)._1;
			if(i < size() - 1) {
				list += ", ";
			}
		}
		return list;
	}

	public double getcRmsd(Tuple2<String, SimplePolymerChain> chain1, Tuple2<String, SimplePolymerChain> chain2, int start1, int end1, int start2, int end2) {
		int s1 = start1;
		int e1 = end1;
		int s2 = start2;
		List<Point3d> coordinates1 = new ArrayList<Point3d>();
		List<Point3d> coordinates2 = new ArrayList<Point3d>();
		for(int n = 0; n < (e1 - s1); n++) {
			Point3d temp1 = chain1._2.getCoordinates()[s1 + n];
			Point3d temp2 = chain2._2.getCoordinates()[s2 + n];
			if((temp1 != null)&&(temp2 != null)) {
				coordinates1.add(temp1);
				coordinates2.add(temp2);
			}
		}
		Point3d[] seq1 = new Point3d[coordinates1.size()];
		Point3d[] seq2 = new Point3d[coordinates2.size()];
		for(int count = 0; count < seq1.length; count ++) {
			seq1[count] = coordinates1.get(count);
			seq2[count] = coordinates2.get(count);
		}
		qcp.set(seq1, seq2);
		return qcp.getRmsd();
	}
}
