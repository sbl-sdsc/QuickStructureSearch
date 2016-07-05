package org.rcsb.projectm.Demo;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;

import scala.Tuple2;

public class Demo1 {

	public static void main(String[] args) {
	// setup spark
			SparkConf conf = new SparkConf()
					.setMaster("local[*]")
					.setAppName("Demo1");	
			
			JavaSparkContext sc = new JavaSparkContext(conf);
			
			String benchmarkDir = "/Users/student/Documents/Data/ProteinResults/EndToEndDistanceSequenceFingerprint_L8B3.7_JaccardIndex_20160627_161208.csv";

			// split input lines of .csv files into chainId1,chainId2, and alignment metrics
			JavaRDD<Float> tm = sc
					.textFile(benchmarkDir, sc.defaultParallelism()) // read files
					.map(t -> t.split(","))
					.map(t -> Float.parseFloat(t[2]))
					.filter(t -> t < 0.5);
			
			long n = tm.count();
			System.out.println(n);
			
			Float sum = tm.reduce((a,b) -> a+b);
			System.out.println("Sum: " + sum);
			System.out.println("Avg: " + sum/n);
			
			JavaRDD<String> map = sc
					.textFile(benchmarkDir, sc.defaultParallelism()) // read files
					.map(t -> t.split(","))
					.mapToPair(t -> new Tuple2<String,String>(t[0],t[1]))
					.map(t -> t._1 + "_" + t._2);
				//	.map(s -> s.split(",")) // split each line into a list items
				//	.cache();
		
			map.foreach(t -> System.out.println(t));
			sc.stop();
			sc.close();
//          ( ͡° ͜ʖ ͡°)
	
	}
}
