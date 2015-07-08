package org.rcsb.project5;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/**
 * This class returns a list of structural clusters given a sequence cluster through an RMSD calculation
 * 
 * @author Justin Li
 * @author Joe Sun
 *
 */
public class StructuralClusterCreator {
	private List<List<Tuple2<String, SimplePolymerChain>>> structuralCluster;
	private SuperPositionQCP qcp;
	private LongestCommonSubstring lcs;

	public StructuralClusterCreator() {
		qcp = new SuperPositionQCP();
		lcs = new LongestCommonSubstring();
	}

	/**
	 * Returns a list of structural clusters given a sequence cluster and a maximum RMSD value for each cluster
	 * @param allClusters
	 * @param maxRmsd
	 * @return structuralCluster
	 */
	public List<List<Tuple2<String, SimplePolymerChain>>> createStructuralCluster(List<Tuple2<String, SimplePolymerChain>> allClusters, double maxRmsd) {
		structuralCluster = new ArrayList<List<Tuple2<String, SimplePolymerChain>>> ();
		for(Tuple2<String, SimplePolymerChain> tempChain: allClusters) {
			boolean unique = true;
			for(List<Tuple2<String, SimplePolymerChain>> cluster: structuralCluster) {
				boolean inCluster = true;
				for(Tuple2<String, SimplePolymerChain> chain: cluster) {
					List<Integer> startEnd = lcs.longestCommonSubstring(tempChain._2.getSequence(), chain._2.getSequence());
					double cRmsd = getcRmsd(tempChain, chain, startEnd.get(0), startEnd.get(1), startEnd.get(2), startEnd.get(3));
//					System.out.println("RMSD of " + tempChain._1 + " and " + chain._1 + ": " + cRmsd);
					if(cRmsd > maxRmsd)
						inCluster = false;
				}
				if(inCluster == true) {
					cluster.add(tempChain);
					unique = false;
				}
			}
			if(unique == true) {
				List<Tuple2<String, SimplePolymerChain>> newCluster = new ArrayList<Tuple2<String, SimplePolymerChain>>();
				newCluster.add(tempChain);
				structuralCluster.add(newCluster);
			}
		}
		return structuralCluster;
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