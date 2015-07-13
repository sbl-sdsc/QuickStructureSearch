package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SQLContext;
import org.apache.spark.sql.SaveMode;
import org.apache.spark.sql.types.DataTypes;
import org.apache.spark.sql.types.StructField;
import org.apache.spark.sql.types.StructType;
import org.rcsb.hadoop.io.HadoopToParquetFile;
import org.rcsb.hadoop.io.HadoopToRowMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.hadoop.io.SimplePolymerChainCodec;
import org.rcsb.hadoop.io.SimplePolymerType;

import scala.Tuple2;

public class HadoopToParquet {

	private static int NUM_THREADS = 8;
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
		
		// read sequence file and map
		JavaRDD<Row> rowRDD = sc
				.sequenceFile(path, Text.class, Text.class)
				//.sample(false, 0.01, 123)
				.mapToPair(t -> new Tuple2<String,String>(t._1.toString(), t._2.toString()))
				.groupByKey()
				.map(new HadoopToParqRow())
				.cache();

		List<StructField> fields = new ArrayList<StructField>();
		fields.add(DataTypes.createStructField("index", DataTypes.StringType, false));
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

		dataFrame.write().mode(SaveMode.Overwrite)
		.partitionBy("index")
		.parquet("/Users/hina/Data/ExampleFiles/seq_whithoutpartion.parquet");
		
		sc.close();
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	private static JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
		.setMaster("local[" + NUM_THREADS + "]")
		.setAppName(HadoopToParquetFile.class.getSimpleName())
		.set("spark.driver.maxResultSize", "2g");
		//.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		//conf.registerKryoClasses(new Class[]{SimplePolymerChain.class, SimplePolymerType.class, SimplePolymerChainCodec.class});

		JavaSparkContext sc = new JavaSparkContext(conf);	
		return sc;
	}
}
