package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.SaveMode;
import org.apache.spark.sql.types.DataTypes;
import org.apache.spark.sql.types.StructField;
import org.apache.spark.sql.types.StructType;

import scala.Tuple2;
/**
 * This class creates a Parquet file from the Hadoop sequence file of Protein-Ligand interactions
 * @author Hinna Shabir
 *
 */
public class HadoopToParquet {

	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
/**
 * 
 * @param args Path of the hadoop sequence file
 * @throws FileNotFoundException
 */
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];
		JavaSparkContext sc = getSparkContext();
		// sc is an existing JavaSparkContext.
		SQLContext sqlContext = new org.apache.spark.sql.SQLContext(sc);
		sqlContext.setConf("spark.sql.parquet.compression.codec", "snappy");
		sqlContext.setConf("spark.sql.parquet.filterPushdown", "true");
		long start = System.nanoTime();
		// read sequence file and map
		JavaRDD<Row> rowRDD = sc
				.sequenceFile(path, Text.class, Text.class)
				//.sample(false, 0.01, 123)
				.mapToPair(t -> new Tuple2<String,String>(t._1.toString(), t._2.toString()))
				.groupByKey()
				.map(new HadoopToParqRow())
				.cache();

		List<StructField> fields = new ArrayList<StructField>();// create data fields of features for the DataFrame 
		fields.add(DataTypes.createStructField("index", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("chainId1", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("chainId2", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("Rnum1", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("Rnum2", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("Ins1", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("Ins2", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("res1", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("res2", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("atom1", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("atom2", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("element1", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("element2", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("distance", DataTypes.IntegerType, false));
		fields.add(DataTypes.createStructField("pdbId", DataTypes.createArrayType(DataTypes.StringType), false));
		StructType schema = DataTypes.createStructType(fields);

		// Apply the schema to the RDD.
		DataFrame dataFrame = sqlContext.createDataFrame(rowRDD, schema);
		dataFrame.coalesce(1).write().mode(SaveMode.Overwrite)
		.partitionBy("index")
		.parquet("/Users/hina/Data/ExampleFiles/seq.parquet");
		sc.close();
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
/**
 * 
 * @return
 */
	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
		.setMaster("local[" + NUM_THREADS + "]")
		.setAppName(HadoopToParquet.class.getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);	
		return sc;
	}
}
