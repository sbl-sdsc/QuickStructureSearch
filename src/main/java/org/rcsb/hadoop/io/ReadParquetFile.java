package org.rcsb.hadoop.io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.Column;
import org.apache.spark.sql.DataFrame;
// Import Row.
import org.apache.spark.sql.Row;
// Import RowFactory.
import org.apache.spark.sql.RowFactory;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.SaveMode;
// Import factory methods provided by DataTypes.
import org.apache.spark.sql.types.DataTypes;
import org.apache.spark.sql.types.StructField;
// Import StructType and StructField
import org.apache.spark.sql.types.StructType;

import scala.collection.mutable.ArrayBuffer;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class ReadParquetFile {    
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main (String[] args) {
		String path = args[0];

		JavaSparkContext sc = getSparkContext();
		// sc is an existing JavaSparkContext.

		SQLContext sqlContext = new org.apache.spark.sql.SQLContext(sc);
		sqlContext.setConf("spark.sql.parquet.filterPushdown", "true");

		long t1 = System.nanoTime();
		DataFrame data = sqlContext.read().parquet(path);

		data.registerTempTable("data");

		DataFrame dataA = sqlContext.sql("SELECT hash, id, chain FROM data where hash IN ('A','B')");
//		DataFrame dataA = sqlContext.sql("SELECT hash, id, chain FROM data where hash IN ('A','B')");

		Row[] rows = dataA.collect();
		for (Row row: rows) {
			ArrayBuffer buffer = row.getAs("chain");
			int[] array = new int[buffer.size()];
			buffer.copyToArray(array);
		//	System.out.println(Arrays.toString(array));
		}
		System.out.println("Records in subset: " + dataA.count());
		
		long t2 = System.nanoTime();
		
		System.out.println("Time: " + (t2-t1)/1E6 + " ms");

	}

	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
		.setMaster("local[" + NUM_THREADS + "]")
		.setAppName(HadoopToParquetFile.class.getSimpleName())
		.set("spark.driver.maxResultSize", "2g")
		.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		return sc;
	}
}
