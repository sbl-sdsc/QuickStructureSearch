package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
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
//		long querystart = System.nanoTime();

		System.out.println("Choose either:  1. Atom name or 2.Element name  3. Simple query");
		Scanner scan = new Scanner(System.in);
		int choice = scan.nextInt();   

		switch (choice) {
		case 3:
			
			System.out.println("Enter the first Protein residue name and atom name: ");
			Scanner s1 = new Scanner(System.in);
			String input = s1.nextLine();   
			String[] str = input.trim().split("\\s+");
			String [] Pro1= new String[2];
			for(int i=0;i<2;i++){
				Pro1[i]="'"+str[i]+"'";
			}
			System.out.println("Enter the first Ligand residue name and atom name: ");
			Scanner s2 = new Scanner(System.in);
			String inp2 = s2.nextLine();   
			String[] st2 = inp2.trim().split("\\s+");
			String [] Lig1= new String[2];
			for(int i=0;i<2;i++){
				Lig1[i]="'"+st2[i]+"'";
			}

			System.out.println("Enter the Distance range  : ");
			Scanner s3 = new Scanner(System.in);
			int [] dist = new int[2];
			for (int i = 0; i < 2; i++) {
				if (s3.hasNextInt()) {
					dist[i]=s3.nextInt();
				}
				else {
					System.out.println("You didn't provide enough numbers");
					System.exit(-1);
				}
			}
			
//			
//			System.out.println("Enter the second Protein residue name and atom name: ");
//			Scanner s11 = new Scanner(System.in);
//			String inp11 = s11.nextLine();   
//			String[] st11 = inp11.trim().split("\\s+");
//			String [] Pro2= new String[2];
//			for(int i=0;i<2;i++){
//				Pro2[i]="'"+st11[i]+"'";
//			}
//			System.out.println("Enter the second Ligand residue name and atom name: ");
//			Scanner s22 = new Scanner(System.in);
//			String inp22 = s22.nextLine();   
//			String[] st22 = inp22.trim().split("\\s+");
//			String [] Lig2= new String[2];
//			for(int i=0;i<2;i++){
//				Lig2[i]="'"+st22[i]+"'";
//			}
//
//			System.out.println("Enter the Distance range  : ");
//			Scanner s33 = new Scanner(System.in);
//			int [] dist2 = new int[2];
//			for (int i = 0; i < 2; i++) {
//				if (s33.hasNextInt()) {
//					dist2[i]=s33.nextInt();
//				}
//				else {
//					System.out.println("You didn't provide enough numbers");
//					System.exit(-1);
//				}
//			}
//			
			
			System.out.println("Case 8");
			results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
					" WHERE D1.res1="+ Pro1[0]+" AND D1.res2=" +Lig1[0]+ 
					" AND D1.atom1="+ Pro1[1]+" AND D1.atom2="+Lig1[1] + " AND D1.distance >="+  dist[0]+ 
					" AND D1.distance <="+ dist[1]);
			
			
//			results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
//					" WHERE D1.res1="+ Pro1[0]+" AND D1.res2=" +Lig1[0]+ 
//					" AND D1.atom1="+ Pro1[1]+" AND D1.atom2="+Lig1[1] + " AND D1.distance >="+  dist[0]+ 
//					" AND D1.distance <="+ dist[1] +" AND D1.pdbId IN" +
//					" (SELECT D2.pdbId FROM Distances D2"+
//					" WHERE D2.res1="+ Pro2[0]+" AND D2.res2=" +Lig2[0]+ 
//					" AND D2.atom1="+ Pro2[1]+" AND D2.atom2="+Lig2[1] + " AND D2.distance >="+  dist2[0]+ 
//					" AND D2.distance <="+ dist2[1]+ ")");
		break;
		
		case 4:
			
/*			System.out.println("Enter the first Protein residue name and element name: ");
			Scanner s11 = new Scanner(System.in);
			String inp = s11.nextLine();   
			String[] st11 = inp.trim().split("\\s+");
			String [] Pro2= new String[2];
			for(int i=0;i<2;i++){
				Pro2[i]="'"+st11[i]+"'";
			}
			System.out.println("Enter the first Ligand residue name and element name: ");
			Scanner s22 = new Scanner(System.in);
			String inp22 = s22.nextLine();   
			String[] st22 = inp22.trim().split("\\s+");
			String [] Lig2= new String[2];
			for(int i=0;i<2;i++){
				Lig2[i]="'"+st22[i]+"'";
			}

			System.out.println("Enter the Distance range  : ");
			Scanner s33 = new Scanner(System.in);
			int [] dist2 = new int[2];
			for (int i = 0; i < 2; i++) {
				if (s33.hasNextInt()) {
					dist2[i]=s33.nextInt();
				}
				else {
					System.out.println("You didn't provide enough numbers");
					System.exit(-1);
				}
			}
			System.out.println("Case 8");
			results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
					" WHERE D1.res1="+ Pro2[0]+" AND D1.res2=" +Lig2[0]+ 
					" AND D1.element1="+ Pro2[1]+" AND D1.element2="+Lig2[1] + " AND D1.distance >="+  dist2[0]+ 
					" AND D1.distance <="+ dist2[1]);
		break;
		*/
		case 1:

			System.out.println("Enter the first Protein residue,atom and number  : ");
			Scanner scan1 = new Scanner(System.in);
			String input1 = scan1.nextLine();   
			String[] str1 = input1.trim().split("\\s+");
			String [] P1= new String[3];
			for(int i=0;i<3;i++){
				P1[i]="'"+str1[i]+"'";
			}
			//System.out.println("P1 :"+ Arrays.toString(P1));
			//System.out.println("choice :"+ P1[2]);
			System.out.println("Enter the first Ligand residue,atom and number  : ");
			Scanner scan2 = new Scanner(System.in);
			String input2 = scan2.nextLine();   
			String[] str2 = input2.trim().split("\\s+");
			String [] L1= new String[3];
			for(int i=0;i<3;i++){
				L1[i]="'"+str2[i]+"'";
			}

			System.out.println("Enter the Distance range  : ");
			Scanner scan5 = new Scanner(System.in);
			int [] input5 = new int[2];
			for (int i = 0; i < 2; i++) {
				if (scan5.hasNextInt()) {
					input5[i]=scan5.nextInt();
				}
				else {
					System.out.println("You didn't provide enough numbers");
					System.exit(-1);
				}
			}
			int [] D1=input5;

			System.out.println("Enter the Second Protein residue,atom and number  : ");
			Scanner scan3 = new Scanner(System.in);
			String input3 = scan3.nextLine();   
			String[] str3 = input3.trim().split("\\s+");
			String [] P2= new String[3];
			for(int i=0;i<3;i++){
				P2[i]="'"+str3[i]+"'";
			}

			System.out.println("Enter the Second Ligand residue,atom and number  : ");
			Scanner scan4 = new Scanner(System.in);
			String input4 = scan4.nextLine();  
			String[] str4 = input4.trim().split("\\s+");
			String [] L2= new String[3];
			for(int i=0;i<3;i++){
				L2[i]="'"+str4[i]+"'";
			}

			System.out.println("Enter the Distance range  : ");
			Scanner scan6 = new Scanner(System.in);
			int [] input6 = new int[2];
			for (int i = 0; i < 2; i++) {
				if (scan6.hasNextInt()) {
					input6[i]=scan6.nextInt();
				}
				else {
					System.out.println("You didn't provide enough numbers");
					System.exit(-1);
				}
			}
			int [] D2=input6;
			// SQL can be run over RDDs that have been registered as tables.
			//1
			if(P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'2'") && L1[2].equalsIgnoreCase("'3'") && L2[2].equalsIgnoreCase("'4'")){
				System.out.println("Case 1");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.index="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1] +" AND D2.index=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="
						+ P2[1]+" AND D2.atom2=" +L2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]
						+" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1) "
						+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2)");
			}
			//2
			else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'2'") && L1[2].equalsIgnoreCase("'3a'") && L2[2].equalsIgnoreCase("'3b'")){
				System.out.println("Case 2");
				long querystart = System.nanoTime();
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
						P2[1]+" AND D2.atom2=" +L2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]
						+ " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"
						+ " AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
				System.out.println("Querying Time: " + (System.nanoTime() - querystart)/1E9 + " sec.");
			}
			//3
			else if (P1[2].equalsIgnoreCase("'1a'") && P2[2].equalsIgnoreCase("'1b'") && L1[2].equalsIgnoreCase("'2'") && L2[2].equalsIgnoreCase("'3'")){
				System.out.println("Case 3");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
						P2[1]+" AND D2.atom2=" +L2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"
						+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2) ");
			}
			//4
			else if (P1[2].equalsIgnoreCase("'1a'") && P2[2].equalsIgnoreCase("'1b'") && L1[2].equalsIgnoreCase("'2a'") && L2[2].equalsIgnoreCase("'2b'")){
				System.out.println("Case 4");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
						P2[1]+" AND D2.atom2=" +L2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			//5
			else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'2'") && L1[2].equalsIgnoreCase("'3'") && L2[2].equalsIgnoreCase("'3'")){
				System.out.println("Case 5");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1= "+ P1[1]+" AND D1.atom2= "+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1]+" AND D2.res1= " + P2[0]+" AND D2.res2= " +L2[0]+ " AND D2.atom1= "+ 
						P2[1]+" AND D2.atom2= " +L2[1] + " AND D2.distance >= "+  D2[0]+ " AND D2.distance <= " + D2[1]+
						" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}	
			//6
			else if (P1[2].equalsIgnoreCase("'1a'") && P2[2].equalsIgnoreCase("'1b'") && L1[2].equalsIgnoreCase("'2'") && L2[2].equalsIgnoreCase("'2'")){
				System.out.println("Case 6");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
						P2[1]+" AND D2.atom2=" +L2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 1a 1b
			}	
			//7
			else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'1'") && L1[2].equalsIgnoreCase("'2'") && L2[2].equalsIgnoreCase("'3'")){
				System.out.println("Case 7");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
						P2[1]+" AND D2.atom2=" +L2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2 )");
			}	
			//8
			else if (P1[2].equalsIgnoreCase("'1'") && P2[2].equalsIgnoreCase("'1'") && L1[2].equalsIgnoreCase("'2a'") && L2[2].equalsIgnoreCase("'2b'")){
				System.out.println("Case 8");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P1[0]+" AND D1.res2=" +L1[0]+ 
						" AND D1.atom1="+ P1[1]+" AND D1.atom2="+L1[1] + " AND D1.distance >="+  D1[0]+ 
						" AND D1.distance <="+ D1[1]+" AND D2.res1=" + P2[0]+" AND D2.res2=" +L2[0]+ " AND D2.atom1="+ 
						P2[1]+" AND D2.atom2=" +L2[1] + " AND D2.distance >="+  D2[0]+ " AND D2.distance <=" + D2[1]+
						" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 2a 2b
			}	
			else{
				System.out.println("NOT FOUND!");
				System.exit(-1);
			}
			break;

		case 2: 

			System.out.println("Enter the first Protein residue,element and number  : ");
			Scanner scan11 = new Scanner(System.in);
			String input11 = scan11.nextLine();   
			String[] str11 = input11.trim().split("\\s+");
			String [] P11= new String[3];
			for(int i=0;i<3;i++){
				P11[i]="'"+str11[i]+"'";
			}
	
			System.out.println("Enter the first Ligand residue,element and number  : ");
			Scanner scan21 = new Scanner(System.in);
			String input21 = scan21.nextLine();   
			String[] str21 = input21.trim().split("\\s+");
			String [] L11= new String[3];
			for(int i=0;i<3;i++){
				L11[i]="'"+str21[i]+"'";
			}

			System.out.println("Enter the Distance range  : ");
			Scanner scan51 = new Scanner(System.in);
			int [] input51 = new int[2];
			for (int i = 0; i < 2; i++) {
				if (scan51.hasNextInt()) {
					input51[i]=scan51.nextInt();
				}
				else {
					System.out.println("You didn't provide enough numbers");
					System.exit(-1);
				}
			}
			int [] D11=input51;

			System.out.println("Enter the Second Protein residue,element and number  : ");
			Scanner scan31 = new Scanner(System.in);
			String input31 = scan31.nextLine();   
			String[] str31 = input31.trim().split("\\s+");
			String [] P21= new String[3];
			for(int i=0;i<3;i++){
				P21[i]="'"+str31[i]+"'";
			}

			System.out.println("Enter the Second Ligand residue,element and number  : ");
			Scanner scan41 = new Scanner(System.in);
			String input41 = scan41.nextLine();  
			String[] str41 = input41.trim().split("\\s+");
			String [] L21= new String[3];
			for(int i=0;i<3;i++){
				L21[i]="'"+str41[i]+"'";
			}

			System.out.println("Enter the Distance range  : ");
			Scanner scan61 = new Scanner(System.in);
			int [] input61 = new int[2];
			for (int i = 0; i < 2; i++) {
				if (scan61.hasNextInt()) {
					input61[i]=scan61.nextInt();
				}
				else {
					System.out.println("You didn't provide enough numbers");
					System.exit(-1);
				}
			}
			int [] D21=input61;
			// SQL can be run over RDDs that have been registered as tables.
			//1
			if(P11[2].equalsIgnoreCase("'1'") && P21[2].equalsIgnoreCase("'2'") && L11[2].equalsIgnoreCase("'3'") && L21[2].equalsIgnoreCase("'4'")){
				System.out.println("Case 2.1");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1="+ P11[1]+" AND D1.element2="+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1] +" AND D2.res1=" + P21[0]+" AND D2.res2=" +L21[0]+ " AND D2.element1="
						+ P21[1]+" AND D2.element2=" +L21[1] + " AND D2.distance >="+  D21[0]+ " AND D2.distance <=" + D21[1]
								+" AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1) "
								+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2)");
			}
			//2
			else if (P11[2].equalsIgnoreCase("'1'") && P21[2].equalsIgnoreCase("'2'") && L11[2].equalsIgnoreCase("'3a'") && L21[2].equalsIgnoreCase("'3b'")){
				System.out.println("Case 2.2");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1="+ P11[1]+" AND D1.element2="+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1]+" AND D2.res1=" + P21[0]+" AND D2.res2=" +L21[0]+ " AND D2.element1="+ 
						P21[1]+" AND D2.element2=" +L21[1] + " AND D2.distance >="+  D21[0]+ " AND D2.distance <=" + D21[1]+
						" AND D1.atom2<>D2.atom2"+ " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			//3
			else if (P11[2].equalsIgnoreCase("'1a'") && P21[2].equalsIgnoreCase("'1b'") && L11[2].equalsIgnoreCase("'2'") && L21[2].equalsIgnoreCase("'3'")){
				System.out.println("Case 2.3");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1="+ P11[1]+" AND D1.element2="+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1]+" AND D2.res1=" + P21[0]+" AND D2.res2=" +L21[0]+ " AND D2.element1="+ 
						P21[1]+" AND D2.element2=" +L21[1] + " AND D2.distance >="+  D21[0]+ " AND D2.distance <=" + D21[1]+
						" AND D1.atom1<>D2.atom1"+ " AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"
						+ " AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2) ");
			}
			//4
			else if (P11[2].equalsIgnoreCase("'1a'") && P21[2].equalsIgnoreCase("'1b'") && L11[2].equalsIgnoreCase("'2a'") && L21[2].equalsIgnoreCase("'2b'")){
				System.out.println("Case 2.4");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1="+ P11[1]+" AND D1.element2="+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1]+" AND D2.res1=" + P21[0]+" AND D2.res2=" +L21[0]+ " AND D2.element1="+ 
						P21[1]+" AND D2.element2=" +L21[1] + " AND D2.distance >="+  D21[0]+ " AND D2.distance <=" + D21[1]+
						" AND D1.atom2<>D2.atom2 AND D1.atom1<>D2.atom1"+ " AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}
			//5
			else if (P11[2].equalsIgnoreCase("'1'") && P21[2].equalsIgnoreCase("'2'") && L11[2].equalsIgnoreCase("'3'") && L21[2].equalsIgnoreCase("'3'")){
				System.out.println("Case 2.5");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1= "+ P11[1]+" AND D1.element2= "+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1]+" AND D2.res1= " + P21[0]+" AND D2.res2= " +L21[0]+ " AND D2.element1= "+ 
						P21[1]+" AND D2.element2= " +L21[1] + " AND D2.distance >= "+  D21[0]+ " AND D2.distance <= " + D21[1]+
						" AND D1.atom2=D2.atom2" + " AND (D1.chainId1<>D2.chainId1 OR D1.Rnum1<>D2.Rnum1)"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");
			}	
			//6
			else if (P11[2].equalsIgnoreCase("'1a'") && P21[2].equalsIgnoreCase("'1b'") && L11[2].equalsIgnoreCase("'2'") && L21[2].equalsIgnoreCase("'2'")){
				System.out.println("Case 2.6");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1="+ P11[1]+" AND D1.element2="+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1]+" AND D2.res1=" + P21[0]+" AND D2.res2=" +L21[0]+ " AND D2.element1="+ 
						P21[1]+" AND D2.element2=" +L21[1] + " AND D2.distance >="+  D21[0]+ " AND D2.distance <=" + D21[1]+
						" AND D1.atom2=D2.atom2 AND D1.atom1<>D2.atom1"+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND D1.chainId2=D2.chainId2 AND D1.Rnum2=D2.Rnum2 AND D1.Ins2=D2.Ins2");// 1a 1b
			}	
			//7
			else if (P11[2].equalsIgnoreCase("'1'") && P21[2].equalsIgnoreCase("'1'") && L11[2].equalsIgnoreCase("'2'") && L21[2].equalsIgnoreCase("'3'")){
				System.out.println("Case 2.7");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1="+ P11[1]+" AND D1.element2="+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1]+" AND D2.res1=" + P21[0]+" AND D2.res2=" +L21[0]+ " AND D2.element1="+ 
						P21[1]+" AND D2.element2=" +L21[1] + " AND D2.distance >="+  D21[0]+ " AND D2.distance <=" + D21[1]+
						" AND D1.atom1=D2.atom1 "+" AND D1.chainId1=D2.chainId1 AND D1.Rnum1=D2.Rnum1 AND D1.Ins1=D2.Ins1"+
						" AND (D1.chainId2<>D2.chainId2 OR D1.Rnum2<>D2.Rnum2 )");
			}	
			//8
			else if (P11[2].equalsIgnoreCase("'1'") && P21[2].equalsIgnoreCase("'1'") && L11[2].equalsIgnoreCase("'2a'") && L21[2].equalsIgnoreCase("'2b'")){
				System.out.println("Case 2.8");
				results = sqlContext.sql ("SELECT D1.pdbId FROM Distances D1"+
						" INNER JOIN Distances D2 ON D2.pdbId=D1.pdbId WHERE D1.res1="+ P11[0]+" AND D1.res2=" +L11[0]+ 
						" AND D1.element1="+ P11[1]+" AND D1.element2="+L11[1] + " AND D1.distance >="+  D11[0]+ 
						" AND D1.distance <="+ D11[1]+" AND D2.res1=" + P21[0]+" AND D2.res2=" +L21[0]+ " AND D2.element1="+ 
						P21[1]+" AND D2.element2=" +L21[1] + " AND D2.distance >="+  D21[0]+ " AND D2.distance <=" + D21[1]+
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

		// The results of SQL queries are DataFrames and support all the normal RDD operations.
		// The columns of a row in the result can be accessed by ordinal.
		List<String> Rows = results.javaRDD().map(new Function<Row, String>() {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			public String call(Row row) {
				return  row.get(0) +" "; //+row.get(1);
				//+" "+row.getString(2)+" "+row.getString(3)+" "+ row.getString(4)+" "
				//+row.getString(5)+" "
				//+ " Distance: "+row.getInt(6) +" PdbIds: "+ row.get(7);
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

	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
		.setMaster("local[" + NUM_THREADS + "]")
		.setAppName(HadoopToParquetFile.class.getSimpleName())
		.set("spark.driver.maxResultSize", "2g");

		JavaSparkContext sc = new JavaSparkContext(conf);	
		return sc;
	}
}
