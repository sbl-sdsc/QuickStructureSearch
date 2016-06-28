package org.rcsb.project10;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;

import scala.Tuple2;

/**
 * This class calculates binary classification scores.
 * @author Peter Rose
 *
 */
public class FingerprintAnalyzer {


	public static void main(String[] args) {
		// setup spark
				SparkConf conf = new SparkConf()
						.setMaster("local[*]")
						.setAppName("FingerprintAnalyzer");	
				
				JavaSparkContext sc = new JavaSparkContext(conf);
				
				String benchmarkDir = args[0];

				// split input lines of .csv files into chainId1,chainId2, and alignment metrics
				JavaPairRDD<String, Double> tmScores = sc
						.textFile(benchmarkDir, sc.defaultParallelism()) // read files as lines
				        .map(t -> t.split(",")) // split String into String[]
				        .mapToPair(t -> new Tuple2<String, Double>(t[0]+"_"+t[1], Double.parseDouble(t[2])));

 //               System.out.println("Benchmark size: " + tmScores.count());
 //               tmScores.foreach(t -> System.out.println(t));
                 
                 String resultsDir = args[1];
              // split input lines of .csv files into chainId1,chainId2, and alignment metrics
 				 JavaPairRDD<String, Double> fpScores = sc
						.textFile(resultsDir, sc.defaultParallelism()) // read files as lines
				        .map(t -> t.split(",")) // split String into String[]
				        .mapToPair(t -> new Tuple2<String, Double>(t[0]+"_"+t[1], Double.parseDouble(t[2])));
                 
 //                System.out.println("Results size: " + fpScores.count());
 //                fpScores.foreach(t -> System.out.println(t));
                 
                 JavaPairRDD<String, Tuple2<Double, Double>> jScores = tmScores.join(fpScores);
                 jScores.foreach(t -> System.out.println(t));
                 
                 // calculate binary classification
                 double threshold = 0.5;
                 JavaPairRDD<String, String> classificationResults = jScores.mapToPair(new BinaryClassificationMapper(threshold)).cache();        
                 classificationResults.foreach(t -> System.out.println(t));
                
                 // calculate binary classification statistics
                 // see: https://en.wikipedia.org/wiki/Precision_and_recall
                 long tp = classificationResults.filter(t -> t._2.equals("TP")).count();
                 long tn = classificationResults.filter(t -> t._2.equals("TN")).count();
                 long fp = classificationResults.filter(t -> t._2.equals("FP")).count();
                 long fn = classificationResults.filter(t -> t._2.equals("FN")).count();
                 
                 sc.stop();
                 sc.close();
 
                 double sensitivity = 1.0 * tp / (double)(tp + fn);
                 double specificity = 1.0 * tn / (double)(tn + fp);
                 double f1Score  = 2.0 * tp / (double)(2*tp + fp + fn);

                 // print binary classification statistics
                 System.out.println("tp: " + tp);
                 System.out.println("tn: " + tn);
                 System.out.println("fp: " + fp);
                 System.out.println("fn: " + fn);
                 System.out.println("Sensitivity: " + sensitivity);
                 System.out.println("Specificity: " + specificity);
                 System.out.println("F1 score: " + f1Score);
	}
}
