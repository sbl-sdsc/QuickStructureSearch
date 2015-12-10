package org.rcsb.hadoop.io;

import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.Comparator;

import javax.vecmath.Point3d;

import me.lemire.integercompression.FastPFOR;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
/* Spark Java programming APIs. It contains the 
 * RDD classes used for Java, as well as the
 * StorageLevels and SparkContext for java.
 */
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.compress.CombinedTransform;
import org.rcsb.compress.DeltaTransform;
import org.rcsb.compress.IntegerTransform;
import org.rcsb.compress.PFORTransform;
import org.rcsb.compress.UnsignedDeltaTransform;
import org.rcsb.structuralSimilarity.GapFilter;

import scala.Tuple2;

/**
 * Demo Map-Reduce program that shows how to read a Hadoop Sequence file and
 * calculate some simple chain statistics
 * @author  Peter Rose
 */
public class HadoopSequenceFileStatistics implements Serializable {    
	private static final long serialVersionUID = 1L;
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3; // Spark recommends 2-3 tasks per thread
	
	public static void main(String[] args ) throws FileNotFoundException
	{
		String path = args[0];

		// This is the default 2 line structure for Spark applications
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(HadoopSequenceFileStatistics.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		
		JavaSparkContext sc = new JavaSparkContext(conf);

		long start = System.nanoTime();
		
//		IntegerTransform transform = new NullOpTransform();
//		IntegerTransform transform = new DeltaTransform();
	//	IntegerTransform transform = new CombinedTransform(new UnsignedDeltaTransform(), new PFORTransform(new FastPFOR()));
		IntegerTransform transform = new UnsignedDeltaTransform();
//		IntegerTransform transform = new DeltaReverseTransform();
//		IntegerTransform transform = new AncientEgyptianDecomposition(new LeGallWavelet());
//		IntegerTransform transform = new CombinedTransform(new NullOpTransform(), new NullOpTransform());
//		IntegerTransform transform = new CombinedTransform(new DeltaTransform(), new AncientEgyptianDecomposition(new LeGallWavelet()));
//		IntegerTransform transform = new PFORTransform(new IntegratedIntCompressor());
		
		// read sequence file and map sequence length to an RDD
		JavaRDD<Integer> len = sc
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.mapToPair(new HadoopToSimpleChainMapperCDF53(transform))
				.filter(t -> t._2.isProtein())
				.map(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
				.filter(new GapFilter(0, 0)) // filter chains with zero gap length and zero gaps
				.map(t -> t._2.length)
				.cache(); // cache since we are using the JavaRDD multiple times below
		
		long chainCount = len.count();
		long minChainLength = len.min(Comparator.naturalOrder());
		long maxChainLength = len.max(Comparator.naturalOrder());
		
		long residueCount = len.reduce((a,b) -> a + b);
		
		sc.stop();
		sc.close();
				
        System.out.println("Total chains         : " + chainCount);
		System.out.println("Total residues       : " + residueCount);	
		System.out.println("Average chain length : " + residueCount/chainCount);
		System.out.println("Minimum chain length : " + minChainLength);
		System.out.println("Maximum chain length : " + maxChainLength);
		
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
}
