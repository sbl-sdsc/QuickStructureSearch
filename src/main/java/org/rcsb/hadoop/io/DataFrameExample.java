package org.rcsb.hadoop.io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
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

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class DataFrameExample {    

	public static void main (String[] args) {

		// This is the default 2 line structure for spark programs in java
		// The spark.executor.memory can only take the maximum java heapspace set by -Xmx
		SparkConf conf = new SparkConf().setMaster("local[1]")
				.setAppName(Demo.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);

		// sc is an existing JavaSparkContext.
		SQLContext sqlContext = new org.apache.spark.sql.SQLContext(sc);

		// Load a text file and convert each line to a JavaBean.
		JavaRDD<String> people = sc.textFile("/Users/peter/Data/ExampleFiles/Fruits.txt");

		// Generate the schema based on the string of schema
		List<StructField> fields = new ArrayList<StructField>();
			fields.add(DataTypes.createStructField("name", DataTypes.StringType, true));
//			fields.add(DataTypes.createStructField("age", DataTypes.StringType, true));
			fields.add(DataTypes.createStructField("age", DataTypes.createArrayType(DataTypes.IntegerType), false));
//		}
		StructType schema = DataTypes.createStructType(fields);

		// Convert records of the RDD (people) to Rows.
		JavaRDD<Row> rowRDD = people.map(
				new Function<String, Row>() {
					private static final long serialVersionUID = 1L;

					public Row call(String record) throws Exception {
						String[] fields = record.split(",");
						System.out.println("Fields: " + Arrays.toString(fields));
						int[] values = new int[fields.length-1];
						for (int i = 1; i < fields.length; i++) {
							values[i-1] = Integer.parseInt(fields[i]);
						}
						return RowFactory.create(fields[0], values);
					}
				});

		// Apply the schema to the RDD.
		DataFrame peopleDataFrame = sqlContext.createDataFrame(rowRDD, schema);

		// Register the DataFrame as a table.
		peopleDataFrame.registerTempTable("people");

		// SQL can be run over RDDs that have been registered as tables.
		DataFrame results = sqlContext.sql("SELECT name,age FROM people");
		
		results.write().mode(SaveMode.Append).save("/Users/peter/Data/ExampleFiles/Fruits.parquet");
//		results.write().json("/Users/peter/Data/ExampleFiles/Fruits.json");
//		results.write().mode(SaveMode.Append).json("/Users/peter/Data/ExampleFiles/Fruits.json");
		DataFrame data = sqlContext.read().format("parquet").load("/Users/peter/Data/ExampleFiles/Fruits.parquet");
		System.out.println("data read;" + data.toString());

		// The results of SQL queries are DataFrames and support all the normal RDD operations.
		// The columns of a row in the result can be accessed by ordinal.
		List<String> names = results.javaRDD().map(new Function<Row, String>() {
			public String call(Row row) {
				return "Name: " + row.getString(0);
			}
		}).collect();
	}
}
