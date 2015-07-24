package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.SaveMode;
import org.rcsb.hadoop.io.HadoopToParquetFile;

public class QueryInteractions {


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
		//results = sqlContext.sql("SELECT * FROM Distances");

		//unique interactions, without distances
		//results = sqlContext.sql("SELECT count(*) FROM Distances GROUP BY res1,res2,atom1,atom2,element1,element2");

		List<String> aminos=Arrays.asList("'ARG'", "'HIS'","'LYS'","'ASP'","'GLU'","'SER'","'THR'","'ASN'","'GLN'",
		"'CYS'","'GLY'","'PRO'","'ALA'","'VAL'","'ILE'","'LEU'","'MET'","'PHE'","'TYR'","'TRP'");

		try {
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
		}

		results.write().mode(SaveMode.Overwrite).save("/Users/hina/Data/ExampleFiles/Distanes_queryResults.parquet");
		System.out.println("Querying Time: " + (System.nanoTime() - querystart)/1E9 + " sec.");

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

//		long count_withdistance= results.distinct().count();
//		System.out.println("unique interactions: " +count_withdistance);
//		long count_withoutdistance= results.count();
//		System.out.println("unique interactions: " +count_withoutdistance);
		sc.close();

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
