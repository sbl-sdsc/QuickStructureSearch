package org.rcsb.REHS;

public class LogSmithWaterman {
	
	public static<T> double distance(int[] seq1, int[] seq2, double open, double extend)
	{
		int m = seq1.length + 1;
		int n = seq2.length + 1;
//		System.out.println(seq1.length +" "+ seq2.length);
		Double[][][][] matrix = new Double[m][n][3][ 11/*Math.min(m-1, n-1)*/]; //guess maximum length of a gap would be 10
		//the first coordinate gives the letter in seq1, the second gives the letter in seq2
		//the third gives whether it is in no gap (0), so previous move was diagonal, if it
		// was coming from up, so a gap in s2, it is (1), and if a gap in s1, it is (2)
		// the fourth coordinate gives how many gaps have been there previous to the current
		
		//thus if the third coordinate is 0, then the fourth coordinate must also be 0
		//initializing the matrix
		// note that also if the third coordinate is not 0, then the fourth coordinate must also
		//not be 0
		
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				for (int k = 0; k < matrix[0][0].length; k++) {
					for (int l = 0; l < matrix[0][0][0].length; l++) {
						matrix[i][j][k][l] = null;
					}
				}
			}
		}
		
		matrix[0][0][0][0] = (double)0;
		for (int i = 1; i < m; i++) {
			matrix[i][0][0][0] = open + extend*Math.log(i);
		}
		
		for (int j = 1; j < n; j++) {
			matrix[0][j][0][0] = open + extend*Math.log(j);
		}
		
		for (int i = 1; i < m; i++) {
			for (int j = 1; j < n; j++) {
				
				// if the ai and bj are to be aligned together - the resulting score
				double min1 = Double.MAX_VALUE;
				for (int k = 0; k < matrix[0][0].length; k++) {
					for (int l = 0; l < matrix[0][0][0].length; l++) {
						if (matrix[i-1][j-1][k][l] != null) {
							min1 = Math.min(min1, matrix[i-1][j-1][k][l]);
						}
					}
				}
				min1 = Math.min(min1, matrix[i-1][j-1][0][0]); // a bit redundant
				min1 += similarity(seq1,seq2,i-1,j-1);
				
				//if b is going to have a gap, so having previous gaps would have
				//third coordinate as 1
				
				double min2 = Double.MAX_VALUE;
				int lvalue2 = 0;
				for (int l = 1; l < matrix[0][0][0].length; l++) {
					if (matrix[i-1][j][1][l] != null && min2 > matrix[i-1][j][1][l] + extend*Math.log((double)(l+1)/l)) {
						min2 = matrix[i-1][j][1][l] + extend*Math.log((double)(l+1)/l);
						lvalue2 = l;
					}
				}
				
				lvalue2 = Math.min(lvalue2, matrix[0][0][0].length-2);
				
				double min3 = Double.MAX_VALUE;
				for (int l = 1; l < matrix[0][0][0].length; l++) {
					if (matrix[i-1][j][2][l] != null) {
						min3 = Math.min(min3, matrix[i-1][j][2][l] + open);
					}
				}
				min3 = Math.min(min3, matrix[i-1][j][0][0] + open);
				
				//if a is going to have a gap
				double min4 = Double.MAX_VALUE;
				for (int l = 1; l < matrix[0][0][0].length; l++) {
					if (matrix[i][j-1][1][l] != null) {
						min4 = Math.min(min4, matrix[i][j-1][1][l] + open);
					}
				}
				min4 = Math.min(min4, matrix[i-1][j][0][0] + open);
				
				double min5 = Double.MAX_VALUE;
				int lvalue5 = 0;
				for (int l = 1; l < matrix[0][0][0].length; l++) {
					if (matrix[i][j-1][2][l] != null && min5 > matrix[i][j-1][2][l] + extend*Math.log((double)(l+1)/l)) {
						min5 = matrix[i][j-1][2][l] + extend*Math.log((double)(l+1)/l);
						lvalue5 = l;
					}
				}
				lvalue5 = Math.min(lvalue5, matrix[0][0][0].length-2);
				
				matrix[i][j][0][0] = min1;
				if (lvalue2 != 0) {
					matrix[i][j][1][lvalue2+1] = min2;
				}
				matrix[i][j][1][1] = min3;
				matrix[i][j][2][1] = min4;
				if (lvalue5 != 0) {
					matrix[i][j][2][lvalue5+1] = min5;
				}
				
//				if (i == j) System.out.println(min1);

				
			}
		}
		
		double value = Double.MAX_VALUE;
		for (int k = 0; k < matrix[0][0].length; k++) {
			for (int l = 0; l < matrix[0][0][0].length; l++) {
				if (matrix[m-1][n-1][k][l] != null) {
					value = Math.min(value, matrix[m-1][n-1][k][l]);
				}
			}
		}
//		System.out.println(matrix[5][5][0][0]);
//		System.out.println(m +" " + n);
		
		double result = value/Math.log(Math.min(seq1.length, seq2.length));
		result = 1- result;
//		System.out.println(seq1.length +" "+seq2.length+" "+value +" " + result);
		
		return 5*result; //way to transform it so that benchmark can judge; should be done better
	}
	
	public static<T> double similarity(int[] seq1, int[] seq2, int i, int j)
	{
		if (seq1[i] == seq2[j] && i < seq1.length && j < seq2.length) return 0;
		else return 1.5;
		
	}
}
