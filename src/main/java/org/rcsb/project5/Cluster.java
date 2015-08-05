package org.rcsb.project5;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralAlignment.SuperPositionQCP;
import org.rcsb.structuralSimilarity.TmScorer;

import scala.Tuple2;
import scala.Tuple3;

/**
 * This class represents a cluster of protein chains with sequence similarity and structural similarity. 
 * It contains information on the sequence cluster Id and structural cluster Id of all of the chains in
 * the cluster and a list of all the protein chains in the structural cluster listed as Tuple2s. It also
 * provides a method for finding the representative chain of the structural cluster.
 * 
 * @author Justin Li
 * @author Joe Sun
 *
 */

public class Cluster {
	private int seqClusterId; // Sequence Cluster Id
	private int strClusterId; // Structural Cluster Id
	private List<Tuple2<String, SimplePolymerChain>> strCluster;
	private Tuple2<String, SimplePolymerChain> repChain; // representative chain
	private List<Tuple3<String, SimplePolymerChain, Double>> strRepChain;
	private SuperPositionQCP qcp;
	private LongestCommonSubstring lcs;
	private double gapPen; // penalty for each gap in the Point3d array
	private double holePen; // penalty for each hole in the Point3d array

	public Cluster(int seqClusterId, int strClusterId, List<Tuple2<String, SimplePolymerChain>> strCluster) {
		this.seqClusterId = seqClusterId;
		this.strClusterId = strClusterId;
		this.strCluster = strCluster;
		strRepChain = new ArrayList<Tuple3<String, SimplePolymerChain, Double>>();
	}

	public Cluster(int seqClusterId, int strClusterId, List<Tuple2<String, SimplePolymerChain>> strCluster, Tuple2<String, SimplePolymerChain> repChain, double gapPen, double holePen) {
		this.seqClusterId = seqClusterId;
		this.strClusterId = strClusterId;
		this.strCluster = strCluster;
		this.repChain = repChain;
		this.gapPen = gapPen;
		this.holePen = holePen;
		qcp = new SuperPositionQCP();
		lcs = new LongestCommonSubstring();
		strRepChain = new ArrayList<Tuple3<String, SimplePolymerChain, Double>>();
	}
	
	public Cluster(int seqClusterId, int strClusterId, List<Tuple2<String, SimplePolymerChain>> strCluster, Tuple2<String, SimplePolymerChain> repChain, double gapPen, double holePen, List<Tuple3<String, SimplePolymerChain, Double>> strRepChains) {
		this.seqClusterId = seqClusterId;
		this.strClusterId = strClusterId;
		this.strCluster = strCluster;
		this.repChain = repChain;
		this.gapPen = gapPen;
		this.holePen = holePen;
		qcp = new SuperPositionQCP();
		lcs = new LongestCommonSubstring();
		strRepChain = strRepChains;
	}

	/**
	 * returns the Sequence Cluster Id of the cluster
	 * @return seqClusterId the sequence cluster Id of the cluster
	 */
	public int getSeqClusterId() {
		return seqClusterId;
	}

	/**
	 * returns the Structural Cluster Id of the cluster
	 * @return strClusterId the structural cluster Id of the cluster
	 */
	public int getStrClusterId() {
		return strClusterId;
	}
	
	/**
	 * sets the sequence cluster Id of the cluster to a given number
	 * @param seqClusterId the sequence cluster Id of the cluster
	 */
	public void setSeqClusterId(int seqClusterId) {
		this.seqClusterId = seqClusterId;
	}
	
	/**
	 * sets the structural cluster Id of the cluster to a given number
	 * @param strClusterId the structural cluster Id of the cluster
	 */
	public void setStrClusterId(int strClusterId) {
		this.strClusterId = strClusterId;
	}

	/**
	 * returns the entire structural cluster as a List of Tuple2s
	 * @return strCluster the structural cluster as a List<Tuple2<String, SimplePolymerChain>>
	 */
	public List<Tuple2<String, SimplePolymerChain>> getStrCluster() {
		return strCluster;
	}

	/**
	 * returns the representative chain of the structural cluster
	 * @return repChain the representative chain of the structural cluster
	 */
	public Tuple2<String, SimplePolymerChain> getRepChain() {
		return repChain;
	}

	public List<Tuple3<String, SimplePolymerChain, Double>> getStrRepChain() {
		return strRepChain;
	}

	/**
	 * returns the number of chains in the structural cluster
	 * @return strCluster.size()
	 */
	public int size() {
		return strCluster.size();
	}

	/**
	 * sets the representative chain of the cluster as the chain given in the argument
	 * @param tuple the desired representative chain of the cluster
	 */
	public void setRepChain(Tuple2<String, SimplePolymerChain> tuple) {
		repChain = tuple;
	}

	/**
	 * adds a structural representative chain to the list of structural representative chains for the cluster
	 * @param tuple the structural representative chain to be added
	 */
	public void addStrRepChain(Tuple3<String, SimplePolymerChain, Double> tuple) {
		strRepChain.add(tuple);
	}

	/**
	 * Finds the representative chain of the structural cluster by using a scoring method which 
	 * adds scores for chain length and subtracts scores for gaps and holes in the chain. The 
	 * chain with the highest score is set as the representative chain. 
	 */
	public void findRepChain() {
		if(isNull()) {
			return;
		}
		double[] scores = new double[size()];
		double tempScore;
		int gaps;
		int holes;
		for(int n = 0; n < size(); n ++) {
			/*			List<Tuple2<String, SimplePolymerChain>> a = getStrCluster();
			Tuple2<String, SimplePolymerChain> b = a.get(n);
			SimplePolymerChain c = b._2;
			if(c == null) {
				System.out.println("SimplePolymerChain null");
				System.out.println(this.getStrCluster().get(n)._1);
				System.exit(-1);
			}
			Point3d[] d = c.getCoordinates();
			double e = d.length;*/
			/*			System.out.println(this);
			System.out.println(this.getStrCluster().get(n)._1);*/
			tempScore = getStrCluster().get(n)._2.getCoordinates().length;
			Point3d[] coordinates = strCluster.get(n)._2.getCoordinates();
			gaps = 0;
			holes = 0;
			boolean newGap = true;
			for(int index = 0; index < coordinates.length; index ++) {
				if(coordinates[index] == null) {
					holes ++;
					if(newGap)
						gaps ++;
					newGap = false;
				} else
					newGap = true;
			}
			tempScore = tempScore - gaps*gapPen - holes*holePen;
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

	/**
	 * Takes a list of chains in the form of Tuple2s and calculates the sums of all the cRMSD values of each
	 * chain with all other respective chains in the structural cluster. The chain with the smallest such sum
	 * is set as the representative chain of the structural cluster.
	 * @param tiebreakerList
	 */
	public void tiebreaker(List<Tuple2<String, SimplePolymerChain>> tiebreakerList) {
		setRepChain(tiebreakerList.get(0));
		/*double[][] array = new double[tiebreakerList.size()][tiebreakerList.size()];
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
		setRepChain(tiebreakerList.get(minIndex));*/
	}

	/**
	 * finds the structural representative chain of a cluster by finding the chain in the cluster which
	 * has the least value of the sums of the RMSD values of that chain and all other chains in the
	 * cluster
	 */
	public void findStrRepChain() {
		double[][] array = new double[strCluster.size()][strCluster.size()];
		List<Integer> startEnd = null;
		double cRmsd;
		for (int outer = 0; outer < array.length - 1; outer++) {
			for (int inner = outer + 1; inner < array.length; inner++) {
				startEnd = lcs.longestCommonSubstring(
						strCluster.get(outer)._2.getSequence(), strCluster.get(inner)._2.getSequence());
				cRmsd = getcRmsd(strCluster.get(outer), strCluster.get(inner), startEnd.get(0),
						startEnd.get(1), startEnd.get(2), startEnd.get(3));
				array[outer][inner] = cRmsd;
			}
		}
		
		double[] simValue = new double[strCluster.size()];
		for (int outer = 0; outer < simValue.length; outer++) {
			for (int inner = 0; inner < simValue.length; inner++) {
				if (inner > outer) {
					if(simValue[outer] < array[outer][inner]) {
						simValue[outer] = array[outer][inner];
					}
				} else if (inner < outer) {
					if(simValue[outer] < array[inner][outer]) {
						simValue[outer] = array[inner][outer];
					}
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
		strRepChain.add(new Tuple3<String, SimplePolymerChain, Double>(strCluster.get(minIndex)._1, strCluster.get(minIndex)._2, simValue[minIndex]));
	}

	public Cluster merge(Cluster c, int seqClusterId, int strClusterId) {
		double[] scores = new double[2];
		double tempScore;
		Point3d[] coordinates;
		int gaps;
		int holes;
		for(int n = 0; n < 2; n ++) {
			if(n == 0) {
				tempScore = this.getRepChain()._2.getCoordinates().length;
				coordinates = this.getRepChain()._2.getCoordinates();
			} else {
				tempScore = c.getRepChain()._2.getCoordinates().length;
				coordinates = c.getRepChain()._2.getCoordinates();
			}
			gaps = 0;
			holes = 0;
			boolean newGap = true;
			for(int index = 0; index < coordinates.length; index ++) {
				if(coordinates[index] == null) {
					holes ++;
					if(newGap)
						gaps ++;
					newGap = false;
				} else
					newGap = true;
			}
			tempScore = tempScore - gaps*gapPen - holes*holePen;
			scores[n] = tempScore;
		}
		Tuple2<String, SimplePolymerChain> newRepChain;
		if(scores[0] >= scores[1]) {
			newRepChain = this.getRepChain();
		} else {
			newRepChain = c.getRepChain();
		}
		strCluster.addAll(c.getStrCluster());
		strRepChain.addAll(c.getStrRepChain());
		return new Cluster(seqClusterId, strClusterId, strCluster, newRepChain, this.gapPen, this.holePen, strRepChain);
	}

	public double similarity() {
		if(isNull()) {
			return -1;
		}else if(size() <= 1) {
			return 0;
		}

		if(getRepChain() == null) {
			findRepChain();
		}

		double sum = 0;
		List<Integer> startEnd = null;

		for(int i = 0; i < size(); i++) {
			startEnd = lcs.longestCommonSubstring(
					getRepChain()._2.getSequence(), strCluster.get(i)._2.getSequence());
			sum += getcRmsd(getRepChain(), strCluster.get(i), startEnd.get(0),
					startEnd.get(1), startEnd.get(2), startEnd.get(3));
		}
		return sum / (size() - 1);
	}

	public double compare(Cluster c) {
		if(isNull() || c.isNull()) {
			return -1;
		}

		double max = 0;
		double temp;

		for(Tuple2<String, SimplePolymerChain> tuple1: strCluster) {
			for(Tuple2<String, SimplePolymerChain> tuple2: c.getStrCluster()) {
				temp = TmScorer.getFatCatTmScore(tuple1._2.getCoordinates(), tuple2._2.getCoordinates())[1];
				if(temp > max) {
					max = temp;
				}
			}
		}

		return max;
	}

	/**
	 * returns a boolean which identifies whether or not the cluster is null
	 * @return true or false
	 */
	public boolean isNull() {
		if(strClusterId == 0) {
			return true;
		}
		return false;
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

	/**
	 * Returns the cRMSD value of a specific portion of the Point3D arrays of the two chains given the starting and ending positions of the specific portions
	 * @param chain1
	 * @param chain2
	 * @param start1
	 * @param end1
	 * @param start2
	 * @param end2
	 * @return qcp.getRmsd()
	 */
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
