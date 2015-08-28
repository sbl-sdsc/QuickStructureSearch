package org.rcsb.ProteinLigandInteractionSearch;
import java.io.FileNotFoundException;
import java.util.List;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;

import scala.Tuple2;
/**
 * Reads the hadoop seqeunce file of protein-ligand interactions
 * @author Hinna Shabir
 *
 */
public class ReadSeqFile { 
	private static int NUM_THREADS = 1;
/**
 * 
 * @param args Path of the sequence file to be read
 * @throws FileNotFoundException
 */
	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = args[0];
		
		long t1 = System.nanoTime();
		ReadSeqFile creator = new ReadSeqFile();
		creator.run(sequenceFileName);
		System.out.println("Time: " + ((System.nanoTime()-t1)/1E9) + " s");
	}
/**
 * 
 * @param sequenceFileName Path of the sequence file to be read
 * @throws FileNotFoundException
 */
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
	
	/**
	 * 
	 * @return
	 */
	private JavaSparkContext getSparkContext() {
		SparkConf conf = new SparkConf()
				.setMaster("local[" + NUM_THREADS + "]")
				.setAppName(this.getClass().getSimpleName());
		JavaSparkContext sc = new JavaSparkContext(conf);
		return sc;
	}
}