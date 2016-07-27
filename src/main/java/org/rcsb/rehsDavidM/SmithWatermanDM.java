package org.rcsb.rehsDavidM;

import java.util.Map;
import java.util.Map.Entry;

/**
 * Class that calculates the Jaccard index between two vectors derived from hashed Features.
 * @author Varkey Alumootil
 *
 */
public class SmithWatermanDM {



	
	/**
	 * Calculates generalized Jaccard distance(index) between two integer sequences.
	 * @param seq1
	 * @param seq2
	 * @param start
	 * @param extend
	 * @return
	 */
	public static<T> double distance(int seq1[], int seq2[], double start, double extend) 
	{
		//seq1 = going down
		//seq2 = going across
		// array to calculate the scores for each route
		// 0 = moving right, 1 = moving down, 2 = moving diagonally
		double [][][] matrix = new double[seq1.length][seq2.length][3];
		for (int i = 0; i < 3 ; i++){
			if (seq1[0] == seq2[0]){
				matrix [0][0][i] = 1;
			}
			else {
				matrix [0][0][i] = 0 - start;
			}
		}
		// calculating square right of the top left
		
		
		if (matrix [0][0][0] != 1){
			matrix[0][1][0] = matrix[0][0][0] - start;
		}
		else{
			matrix [0][1][0] = matrix[0][0][0] - extend;
		}
		
		//calculating square below top left

		if (matrix [0][0][1] != 1){
			matrix[1][0][1] = matrix[0][0][1] - start;
		}
		else{
			matrix [1][0][1] = matrix[0][0][1] - extend;
		}
		
		//calculating diagonal
		if (seq1[1] == seq2[1]){
			matrix[1][1][0] = matrix[1][0][0] - start;
			matrix[1][1][1] = matrix[0][1][1] - start;
			matrix [1][1][2] = matrix[0][0][2] + 1;
		}
		else {
			matrix [1][1][0] = matrix[1][0][1] - start;
			matrix[1][1][1] = matrix [0][1][1] - start;
			matrix [1][1][2] =  matrix [0][0][2] - start;
		}
		//generating 1st row

		for (int i =2; i <seq2.length; i++){
			if (matrix[0][i-1][0] - matrix[0][i-2][0] != 1){
				matrix [0][i][0] = matrix[0][i-1][0] -extend;
			}
			else{
				matrix[0][i][0] = matrix[0][i-1][0] - start;
			}
		}
		//generating 1st column
		for (int i =2; i <seq1.length; i++){
			
			if (matrix[i][0][1] - matrix[i-1][0][1] != 1){
				matrix [i][0][1] = matrix[i-1][0][1] - extend;
			}
			else{
				matrix[i][0][1] = matrix[i-1][0][1] - start;
			}
		}
		//generating 2nd row
		for (int i = 2; i < seq2.length;i++){
			if (seq2[i] == seq1[1]){
				matrix[1][i][2] = Math.max(Math.max(matrix [0][i-1][0], matrix[0][i-1][1]),matrix[0][i-1][2]) + 1;
			}
			else{
				matrix[1][i][2] = Math.max(Math.max(matrix [0][i-1][0], matrix[0][i-1][1]),matrix[0][i-1][2]) - start;
			}
			if (matrix[1][i-1][0] - Math.max(Math.max(matrix [1][i-2][0], matrix[1][i-2][1]),matrix[1][i-2][2])!= 1){
				matrix [1][i][0] = matrix[1][i-1][0] - extend;
			}
			else{
				matrix[1][i][0] =  Math.max(Math.max(matrix [1][i-1][0], matrix[1][i-1][1]),matrix[1][i-1][2]) - start;
			}
			matrix [1][i][1] =  Math.max(Math.max(matrix [0][i][0], matrix[0][i][1]),matrix[0][i][2]) - start;
		}
		//generating 2nd column
		for (int i = 2; i < seq1.length;i++){
			if (seq1[i] == seq2[1]){
				matrix[i][1][2] = Math.max(Math.max(matrix [i-1][0][0], matrix[i-1][0][1]),matrix[i-1][0][2]) + 1;
			}
			else{
				matrix[i][1][2] = Math.max(Math.max(matrix [i-1][0][0], matrix[i-1][0][1]),matrix[i-1][0][2]) - start;
			}
			if (matrix[i-1][1][0] - Math.max(Math.max(matrix [i-2][1][0], matrix[i-2][1][1]),matrix[i-2][1][2])!= 1){
				matrix [i][1][0] = matrix[i-1][1][0] - extend;
			}
			else{
				matrix[i][1][0] =  Math.max(Math.max(matrix [i-1][1][0], matrix[i-1][1][1]),matrix[i-1][1][2]) - start;
			}
			matrix [i][1][1] =  Math.max(Math.max(matrix [i][0][0], matrix[i][0][1]),matrix[i][0][2]) - start;
		}	
		
		//generating the rest
		for (int y = 2; y < seq1.length; y++){
			for (int x = 2; x < seq2.length; x++){
				//diagonal
				if (seq1[y] == seq2[x]){
					matrix[y][x][2] = Math.max(Math.max(matrix [y-1][x-1][0], matrix[y-1][x-1][1]),matrix[y-1][x-1][2]) + 1;
				}
				else if (matrix[y-1][x-1][2] - Math.max(Math.max(matrix [y-2][x-2][0], matrix[y-2][x-2][1]),matrix[y-2][x-2][2])!= 1){
					matrix[y][x][2] = Math.max(Math.max(matrix [y-2][x-2][0], matrix[y-2][x-2][1]),matrix[y-2][x-2][2]) - extend;
				}	
				else{
					matrix[y][x][2] = Math.max(Math.max(matrix [y-1][x-1][0], matrix[y-1][x-1][1]),matrix[y-1][x-1][2]) - start;
				}
				//column
				if (matrix[y-1][x][0] - Math.max(Math.max(matrix [y-2][x][0], matrix[y-2][x][1]),matrix[y-2][x][2])!= 1){
					matrix [y][x][0] = matrix[y-1][x][0] - extend;
				}
				else{
					matrix[y][x][0] =  Math.max(Math.max(matrix [y-1][x][0], matrix[y-1][x][1]),matrix[y-1][x][2]) - start;
				}
				//row
				if (matrix[y][x-1][1] - Math.max(Math.max(matrix [y][x-2][0], matrix[y][x-2][1]),matrix[y][x-2][2])!= 1){
					matrix [y][x][1] = matrix[y][x-1][0] - extend;
				}
				else{
					matrix[y][x][1] =  Math.max(Math.max(matrix [y][x-1][0], matrix[y][x-1][1]),matrix[y][x-1][2]) - start;
				}
			}
		}
		double max =Math.max(Math.max(matrix[seq1.length-1][seq2.length-1][0], matrix[seq1.length-1][seq2.length-1][1]),matrix[seq1.length-1][seq2.length-1][2]);
		double value = 10 * max/(Math.min(seq1.length, seq2.length)- (start + extend *(Math.abs(seq1.length - seq2.length)-1)));
		System.out.println(max + " " + value + " " + seq1.length + " " + seq2.length);
		return value;
	}
		
		
}
