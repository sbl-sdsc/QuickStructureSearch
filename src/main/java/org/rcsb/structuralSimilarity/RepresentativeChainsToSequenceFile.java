package org.rcsb.structuralSimilarity;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.rcsb.utils.BlastClustReader;

import scala.Tuple2;
/**
 * This class creates a Hadoop sequence file for protein chains. By default, it uses the first
 * representative chains from each 40% sequence identity cluster. The Hadoop sequence file
 * uses a delta encoding of the PDB coordinates as well as BZIP2 block level compression.
 * 
 * @author  Peter Rose
 */
public class RepresentativeChainsToSequenceFile {
	private static int NUM_THREADS = 8;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread

	public static void main(String[] args) throws IOException {
		RepresentativeChainsToSequenceFile rep = new RepresentativeChainsToSequenceFile();
		
		long start = System.nanoTime();
		rep.run(args[0], args[1], Integer.parseInt(args[2]));
		
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
	
	public void run(String inputFileName, String outputDirectory, int sequenceIdentity) throws IOException {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String outputFileName = outputDirectory + sequenceIdentity + "_"	+ timeStamp + ".seq";
		
		// setup spark
		SparkConf conf = new SparkConf()
			.setMaster("local[" + NUM_THREADS + "]")
			.setAppName(this.getClass().getSimpleName())
			.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);
		
		// get a map of pdbChainId and cluster id
		BlastClustReader reader = new BlastClustReader(sequenceIdentity);
        final Broadcast<Map<String,Integer>> clusterMap = sc.broadcast(reader.getPdbChainIdClusterMap());
	
        // get chains from file
        JavaPairRDD<String, Point3d[]> chains = sc
        .sequenceFile(inputFileName, Text.class, ArrayWritable.class, NUM_THREADS*NUM_TASKS_PER_THREAD)  // read protein chains
 //       .sample(true,0.1,123456)
        .mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate[]> pairs
        .cache();
        
        long totalChains = chains.count();
        
        // map chains to quality indicator and cluster id, then reduce by quality indicator for each cluster
        JavaRDD<Tuple2<String, Integer>> reps = chains
        .mapToPair(new ChainQualityMapper()) //map to <pdbId.chainId, quality> pairs
        .mapToPair(new ChainToClusterIdMapper(clusterMap)) // map to <clusterId, <pdbId.chainId, quality>> pairs
        .filter(t -> t._1 >= 0) // eliminate any pairs that are not in a cluster
        .reduceByKey(new ReduceByQuality()) // select the best pair for each cluster
        .sortByKey() // sort by chainId
        .values()
        .cache();
       
        // broadcast representative chain ids
        List<String> chainIds = JavaPairRDD.fromJavaRDD(reps)
        		.keys()
        		.collect();
        Broadcast<Set<String>> repbc = sc.broadcast(new HashSet<String>(chainIds));

        // get representative chains
        List<Tuple2<String, Point3d[]>> representatives = chains
        		.filter(new ChainIdFilter<Point3d[]>(repbc))
        		.collect();
        
        long numRepresentatives = representatives.size();      
     
        // write chains to Hadoop sequence file
        PdbToSequenceFileWriter writer = new PdbToSequenceFileWriter(outputFileName);
 //       representatives = representatives.subList(0, 250);
		for (Tuple2<String, Point3d[]> t: representatives) {
			writer.write(t._1, t._2);
		}
	    writer.close();
	    
	    sc.stop();
	    sc.close();
	  
		System.out.println("Total chains: " + totalChains);
		System.out.println("Representative chains: " + numRepresentatives);
	}
}
