package org.rcsb.project10;

import java.util.Arrays;

import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.SequenceFileOutputFormat;
import org.rcsb.mmtf.spark.data.StructureDataRDD;

import scala.Tuple2;

/**
 * This class converts a reduced MMTF Hadoop sequence file 
 * to a Hadoop sequence file of protein chains (as WritableSegments).
 * 
 * @author Peter Rose
 *
 */
public class ReducedStructuresToChains {

	public static void main(String[] args) {
		float maxResolution = 3.0f;
		float maxRfree = 0.3f;
		int minChainLength = 10;
		
		long start = System.nanoTime();

		new StructureDataRDD("/Users/peter/Data/MMTF/reduced")
	    .filterResolution(maxResolution)
		.filterRfree(maxRfree)
		.getJavaRdd()
		.filter(t -> Arrays.asList(t._2.getExperimentalMethods()).contains("X-RAY DIFFRACTION")) // x-ray structures only
		.flatMapToPair(new GappedSegmentGenerator(minChainLength)) // extract protein chains including gaps (missing atoms)
		.mapToPair(t -> new Tuple2<Text, WritableSegment>(new Text(t._1), t._2)) // convert to Text key value
		.saveAsHadoopFile("/Users/peter/Data/MMTF/x-rayChains.seq", Text.class, WritableSegment.class, SequenceFileOutputFormat.class);

		long end = System.nanoTime();
		System.out.println("Time: " + (end-start)/1E9);
	}
}
