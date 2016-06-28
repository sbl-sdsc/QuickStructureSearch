package org.rcsb.projectec;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;

public class Demo1 {

	public static void main(String[] args) {
		
		//System.setProperty("hadoop.home.dir", "C:\\Users\\Emi");
		// setup spark
		SparkConf conf = new SparkConf().setMaster("local[*]").setAppName("Demo1");

		JavaSparkContext sc = new JavaSparkContext(conf);
		String benchmarkDir = "C:\\Users\\Emi\\RCSBstuff\\PairwiseAlignments\\benchmark_20160626_231254.csv\\benchmark_20160626_231254.csv";
		// split input lines of .csv files into chainId1,chainId2, and alignment
		// metrics
		sc.textFile(benchmarkDir, sc.defaultParallelism()).map(s -> s.split(",")).foreach(t -> System.out.print(t)); // split each line into a list items
		// .cache();

	}
}
