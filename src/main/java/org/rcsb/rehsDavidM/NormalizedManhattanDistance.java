package org.rcsb.rehsDavidM;

import java.io.ByteArrayOutputStream;

import java.util.zip.GZIPOutputStream;

/**
 * This class calculates the Normalized Manhattan Distance between two arrays containing the 12 moments describing the protein.
 * 
 * @author  Michael Wang
 * @author David Mao
 */

public class NormalizedManhattanDistance {
	
	public static double distance(double[] s, double[] t) {
		
		double manhattanDis = 0;
		for(int i = 0; i < s.length; i++)
		{
			manhattanDis += (s[i] - t[i]);
		}
		
		return (1.0)/(1.0 + (manhattanDis / 12));
		//return ( 1- (manhattanDis/(Math.pow(s.length,3.25))));
	}
}
