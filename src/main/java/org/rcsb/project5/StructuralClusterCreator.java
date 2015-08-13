package org.rcsb.project5;

import java.util.ArrayList;
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
				//	System.out.println("RMSD of " + tempChain._1 + " and " + chain._1 + ": " + cRmsd);
					if(cRmsd > maxRmsd)
						inCluster = false;
				}
				if(inCluster == true) {
					cluster.add(tempChain);
					unique = false;
					break;
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
	 * Returns a list of structural clusters given a sequence cluster and a maximum RMSD value for each cluster utilizing the triangle inequality
	 * @param allClusters
	 * @param maxRmsd
	 * @return
	 */
	public List<List<Tuple2<String, SimplePolymerChain>>> createStructuralClusterTI(List<Tuple2<String, SimplePolymerChain>> allClusters, double maxRmsd) {
		boolean[][] thresPairs = new boolean[allClusters.size()][allClusters.size()];
		/* Initialization */
		boolean[] done = new boolean[allClusters.size()];
		for(int i = 0; i < done.length; i++) {
			done[i] = false;
		}
		/*for(int i = 0; i < allClusters.size(); i++) {
			if(allClusters.get(i)._1.equals("2IOK.A") || allClusters.get(i)._1.equals("3ERT.A")) {
				System.out.println(allClusters.get(i)._1 + ": " + allClusters.get(i)._2.getSequence());
			}
		}
		System.exit(-1);*/
		/* Iterations */
		for(int i1 = 0; i1 < done.length; i1++) {
			int n_cmp = 0;
			int k = 0;
			int[] index = new int[done.length - i1];
			double[] rmsd = new double[index.length];

			for(int i2 = i1 + 1; i2 < done.length; i2++) {
				if(!done[i2]) {
					index[k] = i2;
					rmsd[k] = RMSD(allClusters.get(i1), allClusters.get(i2));
					n_cmp += 1;
					if(rmsd[k] <= maxRmsd) {
						thresPairs[i1][i2] = true;
						thresPairs[i2][i1] = true;
					}
					k += 1;
				}
			}

			done[i1] = true;
			int n_left = k;
			List<Tuple2<Integer, Double>> index_rmsd = new ArrayList<Tuple2<Integer, Double>>(k);
			for(int count = 0; count < k; count++) {
				index_rmsd.add(new Tuple2<Integer, Double>(index[count], rmsd[count]));
			}

			if(index_rmsd.size() > 0) {
				quickSort(index_rmsd, 0, index_rmsd.size() - 1);
			}

			double ratio = 0;
			for(int j = 0; ((n_cmp == 0 || ((double)(j+1)/n_cmp) > ratio)) && j<n_left; j++) {
				ratio = ((j+1)/n_cmp);

				for(int l=j+1; l<n_left; l++) {
					if ((index_rmsd.get(l)._2 - index_rmsd.get(j)._2) <= maxRmsd) {
						double tempRmsd = RMSD(allClusters.get(index_rmsd.get(l)._1), allClusters.get(index_rmsd.get(j)._1));
						n_cmp += 1;
						if(tempRmsd <= maxRmsd) {
							thresPairs[index_rmsd.get(j)._1][index_rmsd.get(l)._1] = true;
							thresPairs[index_rmsd.get(l)._1][index_rmsd.get(j)._1] = true;
						}
					}
				}
				done[index_rmsd.get(j)._1] = true;
			}
		}
		/*for(int i = 0; i < thresPairs.length; i++) {
			for(int j = 0; j < thresPairs[i].length; j++) {
				if(thresPairs[i][j]) {
					System.out.print("1 ");
				} else {
					System.out.print("0 ");
				}
			}
			System.out.println();
		}*/
		ArrayList<Tuple2<List<Tuple2<String, SimplePolymerChain>>, List<Integer>>> tempStructuralCluster = new ArrayList<Tuple2<List<Tuple2<String, SimplePolymerChain>>, List<Integer>>> ();
		structuralCluster = new ArrayList<List<Tuple2<String, SimplePolymerChain>>> ();
		for(int a = 0; a < allClusters.size(); a++) {
			boolean unique = true;
			for(Tuple2<List<Tuple2<String, SimplePolymerChain>>, List<Integer>> cluster: tempStructuralCluster) {
				boolean inCluster = true;
				for(int b = 0; b < cluster._1.size(); b++) {
					if(!thresPairs[cluster._2.get(b)][a])
						inCluster = false;
				}
				if(inCluster == true) {
					cluster._1.add(allClusters.get(a));
					cluster._2.add(a);
					unique = false;
					break;
				}
			}
			if(unique == true) {
				Tuple2<List<Tuple2<String, SimplePolymerChain>>, List<Integer>> newCluster = new Tuple2<List<Tuple2<String, SimplePolymerChain>>, List<Integer>>(new ArrayList<Tuple2<String, SimplePolymerChain>>(), new ArrayList<Integer>());
				newCluster._1.add(allClusters.get(a));
				newCluster._2.add(a);
				tempStructuralCluster.add(newCluster);
			}
		}
		for(Tuple2<List<Tuple2<String, SimplePolymerChain>>, List<Integer>> cluster: tempStructuralCluster) {
			structuralCluster.add(cluster._1);
		}
		return structuralCluster;
	}

	/**
	 * Returns the cRMSD value of the two chains over the longest common sequence between them
	 * @param tuple
	 * @param tuple2
	 * @return
	 */
	private double RMSD(Tuple2<String, SimplePolymerChain> tuple, Tuple2<String, SimplePolymerChain> tuple2) {
		List<Integer> startEnd = lcs.longestCommonSubstring(tuple._2.getSequence(), tuple2._2.getSequence());
		return getcRmsd(tuple, tuple2, startEnd.get(0), startEnd.get(1), startEnd.get(2), startEnd.get(3));
	}

	void quickSort (List<Tuple2<Integer, Double>> list, int first, int last)
	{
		int g = first, h = last;

		int midIndex = (first + last) / 2;
		double dividingValue = list.get(midIndex)._2;
		do
		{
			while (list.get(g)._2 < dividingValue)
				g++;
			while (list.get(h)._2 > dividingValue)
				h--;
			if (g <= h)
			{
				//swap(list[g], list[h]);
				Tuple2<Integer, Double> temp = list.get(g);
				list.set(g, list.get(h));
				list.set(h, temp);
				g++;
				h--;
			}
		}
		while (g < h);

		if (h > first)
			quickSort (list,first,h);
		if (g < last)
			quickSort (list,g,last);
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