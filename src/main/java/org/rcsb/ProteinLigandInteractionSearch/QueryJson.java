package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
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
	/**
	 * 
	 * @param args path of parquet file and protein-ligand query text file
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args ) throws FileNotFoundException {
		String path = args[0];
		String jsonfile=args[1];   
		List<String> input= ReadJsn(jsonfile);
		String [] Pro = new String[2]; // residue 1 (protein) names
		String [] Lig = new String[2]; // residue 2 (ligand) names
		String [] atom1 = new String[2]; // atom 1 (protein) names
		String [] atom2 = new String[2]; // atom 2 (ligand) names
		String [] elemnt1 = new String[2]; // element 1
		String [] elemnt2 = new String[2]; // element 2
		int [] dist1 = new int [2]; // lower distance bound
		int [] dist2 = new int [2]; // upper distance bound
		String [] pnum = new String [2];
		String [] lnum = new String [2];
		int choice = 0;
		JavaSparkContext sc = getSparkContext();
		// sc is an existing JavaSparkContext.
		SQLContext sqlContext = new org.apache.spark.sql.SQLContext(sc);
		sqlContext.setConf("spark.sql.parquet.compression.codec", "snappy");
		sqlContext.setConf("spark.sql.parquet.filterPushdown", "true");
		long start = System.nanoTime();
		DataFrame data = sqlContext.read().format("parquet").load(path);	
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
		// Register the DataFrame as a table.
		data.registerTempTable("Distances");
		DataFrame results=null;		
		// make data suitable for SQL by enclosing in quotes
		for (int i=0; i<input.size();i++) {
			String [] strng = input.get(i).trim().split("-");
			Pro[i]="'"+strng[0]+"'";
			Lig[i]="'"+strng[1]+"'";
			atom1[i]="'"+strng[2]+"'";
			atom2[i]="'"+strng[3]+"'";
			elemnt1[i]="'"+strng[4]+"'";
			elemnt2[i]="'"+strng[5]+"'";
			dist1[i]=Integer.parseInt(strng[6]);
			dist2[i]=Integer.parseInt(strng[7]);
			// pnum and lnum are used to identify the eight different types within each case
			pnum[i]="'"+strng[8]+"'"; 
			lnum[i]="'"+strng[9]+"'";
			choice=Integer.parseInt(strng[10]);
		}

		switch (choice) {
		case 1: // query using atom name
			if(pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'4'")){
				System.out.println("Case 1");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.index="+ Pro[0]+" AND D1.res2=" +Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0] +" AND D2.index=" + Pro[1]+" AND D2.res2=" + Lig[1] + " AND D2.atom1="
						+ atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]
								+" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1) "
								+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2)");
			}
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3a'") && lnum[1].equalsIgnoreCase("'3b'")){
				System.out.println("Case 2");
				long querystart = System.nanoTime();
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0] + 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1]+ " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]
								+ " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"
								+ " AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
				System.out.println("Querying Time: " + (System.nanoTime() - querystart)/1E9 + " sec.");
			}
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 3");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0] + 
						" AND D1.distance <="+ dist2[0] +" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1]+ " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1] +
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"
						+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2) ");
			}
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
				System.out.println("Case 4");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.atom1="+ atom1[0] +" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0] + 
						" AND D1.distance <="+ dist2[0] +" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1] + " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+  dist1[1] + " AND D2.distance <=" + dist2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 5");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.atom1= "+ atom1[0] +" AND D1.atom2= "+ atom2[0] + " AND D1.distance >="+ dist1[0] + 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1= " + Pro[1]+" AND D2.res2= " +Lig[1]+ " AND D2.atom1= "+ 
						atom1[1]+" AND D2.atom2= " + atom2[1] + " AND D2.distance >= "+ dist1[1] + " AND D2.distance <= " + dist2[1] +
						" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}	
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'2'")){
				System.out.println("Case 6");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1] + " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}	
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 7");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.atom1="+ atom1[0] +" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.atom1="+ 
						atom1[1]+" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1] +
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2 )");
			}	
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
				System.out.println("Case 8");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2="  +Lig[0]+ 
						" AND D1.atom1="+ atom1[0] +" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1]+ " AND D2.atom1="+ 
						atom1[1] +" AND D2.atom2=" + atom2[1] + " AND D2.distance >="+ dist1[1] + " AND D2.distance <=" + dist2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}	
			else {
				System.out.println("NOT FOUND!");
				System.exit(-1);
			}
			break;

		case 2:// query using element name
			if(pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'4'")){
				System.out.println("Case 2.1");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+ elemnt2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0] +" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="
						+ elemnt1[1]+" AND D2.element2=" + elemnt2[1] + " AND D2.distance >="+  dist1[1]+ " AND D2.distance <=" + dist2[1]
								+" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1) "
								+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2)");
			}
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3a'") && lnum[1].equalsIgnoreCase("'3b'")){
				System.out.println("Case 2.2");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+ elemnt2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1]+ " AND D2.element1="+ 
						elemnt1[1]+" AND D2.element2=" + elemnt2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]+
						" AND D1.atom2<>D2.atom2"+ " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 2.3");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+ elemnt2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
						elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]+
						" AND D1.atom1<>D2.atom1"+ " AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"
						+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2) ");
			}
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
				System.out.println("Case 2.4");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+ elemnt2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+ " AND D2.element1="+ 
						elemnt1[1]+" AND D2.element2=" + elemnt2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]+
						" AND D1.atom2<>D2.atom2 AND D1.atom1<>D2.atom1"+ " AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'2'") && lnum[0].equalsIgnoreCase("'3'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 2.5");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.element1= "+ elemnt1[0]+" AND D1.element2= "+ elemnt2[0] + " AND D1.distance >="+  dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1= " + Pro[1]+" AND D2.res2= " + Lig[1]+ " AND D2.element1= "+ 
						elemnt1[1]+" AND D2.element2= " + elemnt2[1] + " AND D2.distance >= "+  dist1[1]+ " AND D2.distance <= " + dist2[1]+
						" AND D1.atom2=D2.atom2" + " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}	
			else if (pnum[0].equalsIgnoreCase("'1a'") && pnum[1].equalsIgnoreCase("'1b'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'2'")){
				System.out.println("Case 2.6");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+ elemnt2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" + Lig[1]+ " AND D2.element1="+ 
						elemnt1[1]+" AND D2.element2=" + elemnt2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]+
						" AND D1.atom2=D2.atom2 AND D1.atom1<>D2.atom1"+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 1a 1b
			}	
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2'") && lnum[1].equalsIgnoreCase("'3'")){
				System.out.println("Case 2.7");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
						" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+ elemnt2[0] +" AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <="+ dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2=" +Lig[1]+" AND D2.element1="+ 
						elemnt1[1]+" AND D2.element2=" +elemnt2[1] + " AND D2.distance >="+ dist1[1]+" AND D2.distance <=" + dist2[1]+
						" AND D1.atom1=D2.atom1 "+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2 )");
			}	
			else if (pnum[0].equalsIgnoreCase("'1'") && pnum[1].equalsIgnoreCase("'1'") && lnum[0].equalsIgnoreCase("'2a'") && lnum[1].equalsIgnoreCase("'2b'")){
				System.out.println("Case 2.8");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ Pro[0]+" AND D1.res2="+ Lig[0]+ 
						" AND D1.element1=" + elemnt1[0]+" AND D1.element2=" + elemnt2[0] + " AND D1.distance >="+ dist1[0]+ 
						" AND D1.distance <=" + dist2[0]+" AND D2.res1=" + Pro[1]+" AND D2.res2="+ Lig[1]+ " AND D2.element1="+ 
						elemnt1[1]+" AND D2.element2=" + elemnt2[1] + " AND D2.distance >="+ dist1[1]+ " AND D2.distance <=" + dist2[1]+
						" AND D1.atom2<>D2.atom2 AND D1.atom1=D2.atom1"+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			else {
				System.out.println("NOT FOUND!");
				System.exit(-1);
			}
			break;

		case 3:
			System.out.println("CASE THREE");
			results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
					" WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
					" AND D1.atom1="+ atom1[0]+" AND D1.atom2="+ atom2[0] + " AND D1.distance >="+ dist1[0]+ 
					" AND D1.distance <="+ dist2[0]);
			break;

		case 4:
			System.out.println("CASE FOUR");
			results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
					" WHERE D1.res1="+ Pro[0]+" AND D1.res2=" + Lig[0]+ 
					" AND D1.element1="+ elemnt1[0]+" AND D1.element2="+ elemnt2[0] + " AND D1.distance >="+ dist1[0]+ 
					" AND D1.distance <="+ dist2[0]);
			break;

		default:
			System.out.println(" Wrong choice!");
		}
		
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
	 * @return queries string containing parameters for the queries
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
		int choice=0;
		List<String> queries = new ArrayList <String>();
		JSONParser parser = new JSONParser();
		try {
			Object obj = parser.parse(new FileReader(path));
			JSONObject jsonObject = (JSONObject) obj;
			JSONArray interactions = (JSONArray) jsonObject.get("interactions");
			Iterator<?> iterator = interactions.iterator();
			int count= interactions.size();
			while (iterator.hasNext()) {
				JSONObject innerobj = (JSONObject) iterator.next();
				int num1, num2;
				int site1 = Integer.parseInt((String) innerobj.get("res1"));
				int site2 = Integer.parseInt((String) innerobj.get("res2"));
				distmin =Integer.parseInt((String) innerobj.get("distanceMin"));
				distmax =Integer.parseInt((String) innerobj.get("distanceMax"));
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
				if (count==1) {	// Simple query, i.e just one protein-ligand pair
					if (atm1==null && atm2==null)
						choice = 4;		// element name is provided in query
					else if (elmnt1==null && elmnt2==null)
						choice = 3;		// atom name is provided in query
				}
				else if (atm1==null && atm2==null) {
					choice = 2; // i.e. element name is provided in query
				}
				else if (elmnt1==null && elmnt2==null) {
					choice = 1; // atom name is provided in query
				}
				else {
					choice = 0;
				}
				str= res1+"-"+res2+"-"+atm1+"-"+atm2+"-"+elmnt1+"-"+elmnt2+"-"+distmin+"-"+distmax+'-'+num1+'-'+num2+'-'+ choice;	
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