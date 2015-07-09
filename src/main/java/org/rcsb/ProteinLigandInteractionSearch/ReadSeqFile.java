package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileNotFoundException;
import java.util.List;

import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;

import scala.Tuple2;

/**
 * This class creates structural alignments between random protein chain pairs 
 * using jFatCAT and scores the alignments with the TM score
 * 
 * @author  Peter Rose
 */
public class ReadSeqFile { 
	private static int NUM_THREADS = 1;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0]; 
		
		long t1 = System.nanoTime();
		ReadSeqFile creator = new ReadSeqFile();
		creator.run(sequenceFileName);
		System.out.println("Time: " + ((System.nanoTime()-t1)/1E9) + " s");
	}

	private void run(String sequenceFileName) throws FileNotFoundException {
		// setup spark
		JavaSparkContext sc = getSparkContext();
		
		// Get <Interaction, pdbId> pairs
        List<Tuple2<String, String>> interactions = sc
				.sequenceFile(sequenceFileName, Text.class, Text.class, NUM_THREADS)  // read protein chains
				.map(t -> new Tuple2<String,String>(t._1.toString(), t._2.toString()))
				.collect(); // return results to master node
	   for(Tuple2<String, String> interaction: interactions){
		   System.out.println("Interaction: "+ interaction._1 + "  "+ interaction._2);
	   }

		sc.stop();
		sc.close();

	}
	
	private JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName())
				.set("spark.driver.maxResultSize", "2g")
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		return sc;
	}

}

