package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SQLContext;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.rcsb.hadoop.io.HadoopToParquetFile;
/**
 * 
 * @author Hinna Shabir
 *
 */
public class QueryJson {

	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	@SuppressWarnings("null")
	/**
	 * 
	 * @param args
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];
		String jsonfile=args[1];   
		List<String> input= ReadJsn(jsonfile);
		String [] Pro = new String[2]; // residue 1 (protein) names
		String [] Lig= new String[2]; // residue 2 (ligand) names
		String [] atom1= new String[2]; // atom 1 (protein) names
		String [] atom2= new String[2]; // atom 2 (ligand) names
		String [] elemnt1= new String[2]; // element 1
		String [] elemnt2= new String[2]; // element 2
		int [] dist1 = new int [2]; // lower bound
		int [] dist2= new int [2];  // upper bound
		String [] pnum= new String [2];
		String [] lnum= new String [2];
		JavaSparkContext sc = getSparkContext();
		// sc is an existing JavaSparkContext.
		SQLContext sqlContext = new org.apache.spark.sql.SQLContext(sc);
		sqlContext.setConf("spark.sql.parquet.compression.codec", "snappy");
		sqlContext.setConf("spark.sql.parquet.filterPushdown", "true");
		long start = System.nanoTime();

		DataFrame data = sqlContext.read().format("parquet").load(path);
		System.out.println("data read;" + data.toString());
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
		// Register the DataFrame as a table.
		data.registerTempTable("Distances");
		DataFrame results=null;
		System.out.println("Choose either:  1. Atom name or 2.Element name  3. Simple query");
		Scanner scan = new Scanner(System.in);
		int choice = scan.nextInt();
		// make data suitable for sql by enclosing in quotes
	    for (int i=0; i<input.size();i++)
		{
		String [] strng = input.get(i).trim().split("-");
		Pro[i]="'"+strng[0]+"'";
		Lig[i]="'"+strng[1]+"'";
		atom1[i]="'"+strng[2]+"'";
		atom2[i]="'"+strng[3]+"'";
		elemnt1[i]="'"+strng[4]+"'";
		elemnt2[i]="'"+strng[5]+"'";
		dist1[i]=Integer.parseInt(strng[6]);
		dist2[i]=Integer.parseInt(strng[7]);
		pnum[i]="'"+strng[8]+"'";
		lnum[i]="'"+strng[9]+"'";
		}
		switch (choice) {
		case 3:
//			String [] strng=input.get(1).trim().split("-");;
//			Pro[0]="'"+strng[0]+"'";
//			Lig[0]="'"+strng[1]+"'";
//			atom1[0]="'"+strng[2]+"'";
//			atom2[0]="'"+strng[3]+"'";
//			dist1[0]=Integer.parseInt(strng[4]);
//			dist2[0]=Integer.parseInt(strng[5]);
			results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
					" WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
					" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+  dist1[0]+ 
					" AND D1.distance <="+ dist2[0]);
			break;

			case 1:

			if(pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'4'")){
				System.out.println("Case 1");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.index="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+atom2[0] + " AND D1.distance >="+  dist1[0]+ 
						" AND D1.distance <="+ dist2[0] +" AND D2.index=" + Pro[1]+" AND D2.res2=" + Lig[1] + " AND D2.atom1="
						+ atom1[1]+" AND D2.atom2=" +atom2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]
						+" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1) "
						+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2)");
			}
			//2
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3a'") && lnum[1].equalsIgnoreCase("'3b'")){
				System.out.println("Case 2");
				long querystart = System.nanoTime();
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+atom2[0] + " AND D1.distance >="+ dist1[0] + 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]
						+ " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"
						+ " AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
				System.out.println("Querying Time: " + (System.nanoTime() - querystart)/1E9 + " sec.");
			}
			//3
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 3");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+atom2[0] + " AND D1.distance >="+ dist1[0] + 
						" AND D1.distance <="+ dist2[0] +" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1]+ " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1] +
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"
						+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2) ");
			}
			//4
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
				System.out.println("Case 4");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0] +" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+  dist1[0] + 
						" AND D1.distance <="+ dist2[0] +" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1] + " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+  dist1[1] + " AND D2.distance <=" + dist2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			//5
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'3'")){
				long Qstart = System.nanoTime();
				System.out.println("Case 5");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1= "+ atom1[0] +" AND D1.atom2= "+ atom2[0] + " AND D1.distance >="+ dist1[0] + 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1= " + Pro[1]+" AND D2.res2= " +Lig[1]+ " AND D2.atom1= "+ 
						atom1[1]+" AND D2.atom2= " + atom2[1] + " AND D2.distance >= "+ dist1[1] + " AND D2.distance <= " + dist2[1] +
						" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
				System.out.println("Querying Time: " + (System.nanoTime() - Qstart)/1E9 + " sec.");
			}	
			//6
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'2'")){
				System.out.println("Case 6");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+  dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1] + " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 1a 1b
			}	
			//7
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 7");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0] +" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+  dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1] +
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2 )");
			}	
			//8
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
				System.out.println("Case 8");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0] +" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+  dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.atom1="+ 
						atom1[1] +" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+  dist1[1] + " AND D2.distance <=" + dist2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 2a 2b
			}	
			else{
				System.out.println("NOT FOUND!");
				System.exit(-1);
			}
			break;
			case 2:
				if(pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'4'")){
					System.out.println("Case 2.1");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0] +" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="
							+ elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]
									+" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1) "
									+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2)");
				}
				//2
				else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3a'") && lnum[1].equalsIgnoreCase("'3b'")){
					System.out.println("Case 2.2");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
							elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
							" AND D1.atom2<>D2.atom2"+ " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
							" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
				}
				//3
				else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
					System.out.println("Case 2.3");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
							elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
							" AND D1.atom1<>D2.atom1"+ " AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"
							+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2) ");
				}
				//4
				else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
					System.out.println("Case 2.4");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
							elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
							" AND D1.atom2<>D2.atom2 AND D1.atom1<>D2.atom1"+ " AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
							" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
				}
				//5
				else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'3'")){
					System.out.println("Case 2.5");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1= "+ elemnt1[0]+" AND D1.element2= "+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0]+" AND D2.res1= " + Pro[1]+" AND D2.res2= " +Lig[1]+ " AND D2.element1= "+ 
							elemnt1[1]+" AND D2.element2= " +elemnt2[1] + " AND D2.distance >= "+  dist1[1]+ " AND D2.distance <= " + dist2[1]+
							" AND D1.atom2=D2.atom2" + " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
							" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
				}	
				//6
				else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'2'")){
					System.out.println("Case 2.6");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
							elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
							" AND D1.atom2=D2.atom2 AND D1.atom1<>D2.atom1"+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
							" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 1a 1b
				}	
				//7
				else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
					System.out.println("Case 2.7");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
							elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
							" AND D1.atom1=D2.atom1 "+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
							" AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2 )");
				}	
				//8
				else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
					System.out.println("Case 2.8");
					results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
							" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
							" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
							" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
							elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
							" AND D1.atom2<>D2.atom2 AND D1.atom1=D2.atom1"+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
							" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 2a 2b
				}

				else{
					System.out.println("NOT FOUND!");
					System.exit(-1);
				}

				break;
		default:
			System.out.println("Choose either 1 or 2");
		}
//		System.out.println("Querying Time: " + (System.nanoTime() - querystart)/1E9 + " sec.");
		List<String> Rows = results.javaRDD().map(new Function<Row, String>() {
			private static final long serialVersionUID = 1L;
			public String call(Row row) {
				return  row.get(0) +" "; 
			}
		}).collect();
		if (Rows.isEmpty()){
			System.out.println("QUERY RETURNS NO RESULTS");
		}
		for (String s: Rows){
			System.out.println(s.toCharArray());
		}

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
/**
 * 
 * @param path
 * @return
 */
	private static List<String> ReadJsn (String path) {
		String res1=null;
		String res2= null;
		String atm1=null;
		String atm2=null;
		String elmnt1=null;
		String elmnt2=null;
		int distmin = 0;
		int distmax = 0;
		String str=null;
		List<String> queries = new ArrayList <String>();
		JSONParser parser = new JSONParser();
		try {
//			String path=args[0];
			Object obj = parser.parse(new FileReader(path));
			JSONObject jsonObject = (JSONObject) obj;
			JSONArray interactions = (JSONArray) jsonObject.get("interactions");
			Iterator<?> iterator = interactions.iterator();
			while (iterator.hasNext()) {
				JSONObject innerobj = (JSONObject) iterator.next();
//				System.out.println("res1 :"+innerobj.get("res1")+"res2 :"+innerobj.get("res2")+"distance:"+innerobj.get("distanceMin")+"-"+innerobj.get("distanceMax"));
				int site1 = Integer.parseInt((String) innerobj.get("res1"));
				int site2 = Integer.parseInt((String) innerobj.get("res2"));
				distmin =Integer.parseInt((String) innerobj.get("distanceMin"));
				distmax =Integer.parseInt((String) innerobj.get("distanceMax"));
				int num1, num2;
//				System.out.println("1: "+ site1);System.out.println("2: "+ site2);
				JSONArray sites= (JSONArray) jsonObject.get("sites");
				JSONObject innerObject1 =(JSONObject) sites.get(site1-1);
				JSONObject innerObject2 =(JSONObject) sites.get(site2-1);
				res1 =  (String) innerObject1.get("residueName");
				res2 =  (String) innerObject2.get("residueName");
				atm1 = (String) innerObject1.get("atomName");
				atm2 = (String) innerObject2.get("atomName");
				elmnt1 = (String) innerObject1.get("element");
				elmnt2 = (String) innerObject2.get("element");
				num1= Integer.parseInt((String) innerObject1.get("num"));
				num2= Integer.parseInt((String) innerObject2.get("num"));
				System.out.println(atm1);
				System.out.println(atm2);
				str= res1+"-"+res2+"-"+atm1+"-"+atm2+"-"+elmnt1+"-"+elmnt2+"-"+distmin+"-"+distmax+'-'+num1+'-'+num2;	
				System.out.println("String: "+str);
				queries.add(str);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		return (queries);
	}
	/**
	 * 
	 * @return
	 */

	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
		.setMaster("local[" + NUM_THREADS + "]")
		.setAppName(HadoopToParquetFile.class.getSimpleName())
		.set("spark.driver.maxResultSize", "2g");

		JavaSparkContext sc = new JavaSparkContext(conf);	
		return sc;
	}
}
