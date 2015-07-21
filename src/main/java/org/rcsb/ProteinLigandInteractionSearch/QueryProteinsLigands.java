package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SQLContext;
import org.rcsb.hadoop.io.HadoopToParquetFile;

public class QueryProteinsLigands {


	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];
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
		long querystart = System.nanoTime();
		// SQL can be run over RDDs that have been registered as tables.

		//results = sqlContext.sql("SELECT count(*) from Distances");
		//unique interactions, with distances.
		//	results = sqlContext.sql("SELECT * FROM Distances");

		//unique interactions, without distances
		//results = sqlContext.sql("SELECT count(*) FROM Distances GROUP BY res1,res2,atom1,atom2,element1,element2");

		//List<String> aminos=Arrays.asList("'ARG'", "'HIS'","'LYS'","'ASP'","'GLU'","'SER'","'THR'","'ASN'","'GLN'",
		//"'CYS'","'GLY'","'PRO'","'ALA'","'VAL'","'ILE'","'LEU'","'MET'","'PHE'","'TYR'","'TRP'");

		System.out.println("Enter the first Protein residue,atom and number  : ");
		Scanner scan = new Scanner(System.in);
		String input1 = scan.nextLine();   
		String[] str1 = input1.trim().split("\\s+");
		String [] P1= new String[3];
		for(int i=0;i<3;i++){
			P1[i]="'"+str1[i]+"'";
		}
		//System.out.println("P1 :"+ Arrays.toString(P1));
		//System.out.println("choice :"+ P1[2]);
		System.out.println("Enter the first Ligand residue,atom and number  : ");
		Scanner scan2 = new Scanner(System.in);
		String input2 = scan.nextLine();   
		String[] str2 = input2.trim().split("\\s+");
		String [] L1= new String[3];
		for(int i=0;i<3;i++){
			L1[i]="'"+str2[i]+"'";
		}
		
		System.out.println("Enter the Distance range  : ");
		Scanner scan5 = new Scanner(System.in);
		int [] D1 = new int[2];
		for (int i = 0; i < 2; i++) {
		     if (scan5.hasNextInt()) {
		        D1[i]=scan5.nextInt();
		     }
		     else {
	                System.out.println("You didn't provide enough numbers");
	                break;
	            }
		  }
		
	/*	for(int i=0;i<2;i++){
			D1[i]="'"+D1[i]+"'";
		}*/
		
		System.out.println("Enter the Second Protein residue,atom and number  : ");
		Scanner scan3 = new Scanner(System.in);
		String input3 = scan.nextLine();   
		String[] str3 = input3.trim().split("\\s+");
		String [] P2= new String[3];
		for(int i=0;i<3;i++){
			P2[i]="'"+str3[i]+"'";
		}
		
		System.out.println("Enter the Second Ligand residue,atom and number  : ");
		Scanner scan4 = new Scanner(System.in);
		String input4 = null;
		//if (scan4.hasNextLine()) {
		input4 = scan.nextLine();  
		//}
		//else{
		//	System.out.println("You didn't provide enough input");
		//}
		String[] str4 = input4.trim().split("\\s+");
		String [] L2= new String[3];
		for(int i=0;i<3;i++){
			L2[i]="'"+str4[i]+"'";
		}
		
		System.out.println("Enter the Distance range  : ");
		Scanner scan6 = new Scanner(System.in);
		int [] D2 = new int[2];
		for (int i = 0; i < 2; i++) {
		     if (scan6.hasNextInt()) {
		        D2[i]=scan6.nextInt();
		     }
		     else {
	                System.out.println("You didn't provide enough numbers");
	                break;
	            }
		  }
		
		//1
		if(P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'2'") && L1[2].equalsIgnoreCase("'3'") && L2[2].equalsIgnoreCase("'4'")){
			System.out.println("Case 1");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					"AND D1.chainId1!=D2.chainId1 AND D1.Rnum1!=D2.Rnum1 AND D1.Ins1!=D2.Ins1 AND "
					+ "D1.chainId2!=D2.chainId2 AND D1.Rnum2!=D2.Rnum2 AND D1.Ins2!=D2.Ins2");
		}
		//2
		else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'2'") && L1[2].equalsIgnoreCase("'3a'") && L2[2].equalsIgnoreCase("'3b'")){
			System.out.println("Case 2");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					"AND D1.chainId1!=D2.chainId1 AND D1.Rnum1!=D2.Rnum1 AND D1.Ins1!=D2.Ins1"
					+ "AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
		}
		//3
		else if (P1[2].equalsIgnoreCase("'1a'") && P2[2].equalsIgnoreCase("'1b'") && L1[2].equalsIgnoreCase("'2'") && L2[2].equalsIgnoreCase("'3'")){
			System.out.println("Case 3");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					"AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"
					+ "AND D1.chainId2!=D2.chainId2 AND D1.Rnum2!=D2.Rnum2 AND D1.Ins2!=D2.Ins2");
		}
		//4
		else if (P1[2].equalsIgnoreCase("'1a'") && P2[2].equalsIgnoreCase("'1b'") && L1[2].equalsIgnoreCase("'2a'") && L2[2].equalsIgnoreCase("'2b'")){
			System.out.println("Case 4");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
					" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
		}
		//5
		else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'2'") && L1[2].equalsIgnoreCase("'3'") && L2[2].equalsIgnoreCase("'3'")){
			System.out.println("Case 5");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					" AND D1.chainId1!=D2.chainId1 AND D1.Rnum1!=D2.Rnum1 AND D1.Ins1!=D2.Ins1"+
					" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
		}	
		//6
		else if (P1[2].equalsIgnoreCase("'1a'") && P2[2].equalsIgnoreCase("'1b'") && L1[2].equalsIgnoreCase("'2'") && L2[2].equalsIgnoreCase("'2'")){
			System.out.println("Case 6");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
					" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 1a 1b
		}	
		//7
		else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'1'") && L1[2].equalsIgnoreCase("'2'") && L2[2].equalsIgnoreCase("'3'")){
			System.out.println("Case 7");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
					" AND D1.chainId2!=D2.chainId2 AND D1.Rnum2!=D2.Rnum2 AND D1.Ins2!=D2.Ins2");
		}	
		//8
		else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'1'") && L1[2].equalsIgnoreCase("'2a'") && L2[2].equalsIgnoreCase("'2b'")){
			System.out.println("Case 8");
			results = sqlContext.sql ("SELECT D2.pdbId FROM Distances D1"+
					" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
					" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
					" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
					P2[1]+" AND D2.atom2=" +P2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
					" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
					" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 2a 2b
		}
		else{
			System.out.println("NOT FOUND!");
			System.exit(-1);
		}
		

		/*try {
			FileWriter writer = new FileWriter("/Users/hina/Data/#interactions.csv");

		for (int i=0; i<20;i++){
			String protein=aminos.get(i);
			results = sqlContext.sql("SELECT * FROM Distances WHERE res1="+protein);
			long numberofint= results.distinct().count();
			System.out.println(protein+" " +numberofint);
			String count= " "+numberofint;
			writer.append(protein);
			writer.append(",");
		    writer.append(count);
			writer.append("\n");
		}
		writer.flush();
	    writer.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}*/
		//'ASP' '017' 'OD2' 'O18' '24' '26'
		//'GLY' '017' 'O' 'N20' '30' '32'
		//results.write().mode(SaveMode.Overwrite).save("/Users/hina/Data/ExampleFiles/Distanes_queryResults.parquet");
		System.out.println("Querying Time: " + (System.nanoTime() - querystart)/1E9 + " sec.");

		// The results of SQL queries are DataFrames and support all the normal RDD operations.
		// The columns of a row in the result can be accessed by ordinal.
		List<String> Rows = results.javaRDD().map(new Function<Row, String>() {
			public String call(Row row) {
				return  row.get(0) +" "; //+row.get(1);
				//+" "+row.getString(2)+" "+row.getString(3)+" "+ row.getString(4)+" "
				//+row.getString(5)+" "
				//+ " Distance: "+row.getInt(6) +" PdbIds: "+ row.get(7);
			}
		}).collect();

		for (String s: Rows){
			System.out.println(s);
		}
		/*		long count_withdistance= results.distinct().count();
		System.out.println("unique interactins: " +count_withdistance);
		long count_withoutdistance= results.count();
		System.out.println("unique interactins: " +count_withoutdistance);
		sc.close();*/

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
		.setMaster("local[" + NUM_THREADS + "]")
		.setAppName(HadoopToParquetFile.class.getSimpleName())
		.set("spark.driver.maxResultSize", "2g");

		JavaSparkContext sc = new JavaSparkContext(conf);	
		return sc;
	}
}
