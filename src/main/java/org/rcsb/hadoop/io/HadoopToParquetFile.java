package org.rcsb.hadoop.io;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.SaveMode;
import org.apache.spark.sql.types.DataTypes;
import org.apache.spark.sql.types.StructField;
import org.apache.spark.sql.types.StructType;

import scala.Tuple2;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class HadoopToParquetFile {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];

		JavaSparkContext sc = getSparkContext();
		// sc is an existing JavaSparkContext.
		
		SQLContext sqlContext = new org.apache.spark.sql.SQLContext(sc);
//		sqlContext.setConf("spark.sql.parquet.compression.codec", "gzip");
//		sqlContext.setConf("spark.sql.parquet.compression.codec", "uncompressed");
		sqlContext.setConf("spark.sql.parquet.compression.codec", "snappy");
		sqlContext.setConf("spark.sql.parquet.filterPushdown", "true");
		long start = System.nanoTime();
		
		// if you need both the coordinates and the sequences, use this section of code
		// read sequence file and map to PdbId.chainId, SimplePolymerChain pairs
		JavaRDD<Row> rowRDD = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
//				.sample(false, 0.01, 123)
				.map(new HadoopToRowMapper())
				.cache();

		List<StructField> fields = new ArrayList<StructField>();
		fields.add(DataTypes.createStructField("hash", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("id", DataTypes.StringType, false));
		fields.add(DataTypes.createStructField("chain", DataTypes.createArrayType(DataTypes.IntegerType), false));
		StructType schema = DataTypes.createStructType(fields);
		
		// Apply the schema to the RDD.
		DataFrame dataFrame = sqlContext.createDataFrame(rowRDD, schema);
//
//		// Register the DataFrame as a table.
//		dataFrame.registerTempTable("chains");
//		
//		int n = 0;
//		long t1 = System.nanoTime();
//		List<Row> sample = rowRDD.sample(false,  0.01, 123).collect();
//		StringBuilder sb = new StringBuilder();
//		sb.append("(");
//        for (Row r: sample) {
//        	sb.append("'");
//        	sb.append(r.getString(0));
//        	sb.append("',");
// //       	System.out.println("sample: " + r.getString(0));
//        	n++;
//        };
//        sb.setCharAt(sb.length()-1, ')');
//        System.out.println(sb.toString());
//        
//		// SQL can be run over RDDs that have been registered as tables.
//	    DataFrame results = sqlContext.sql("SELECT id FROM chains where id IN " + sb.toString());
//
//        long t2 = System.nanoTime();
//        System.out.println("Looked up: " + n);
//        System.out.println("Time: " + (t2-t1)/1E6 + " ms");
//        System.out.println("Time per entry: " + (t2-t1)/1E6/n + " ms");
		
		dataFrame.write().mode(SaveMode.Overwrite).partitionBy("hash").parquet("/Users/peter/Data/ExampleFiles/chains.parquet");

		sc.close();

		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(HadoopToParquetFile.class.getSimpleName())
				.set("spark.driver.maxResultSize", "2g")
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		conf.registerKryoClasses(new Class[]{SimplePolymerChain.class, SimplePolymerType.class, SimplePolymerChainCodec.class});

		JavaSparkContext sc = new JavaSparkContext(conf);
//		sc.hadoopConfiguration().setInt("parquet.block.size", 1012*1024*256);
//		sc.hadoopConfiguration().setInt("dfs.block.size", 1012*1024*256);
		
		return sc;
	}
}
