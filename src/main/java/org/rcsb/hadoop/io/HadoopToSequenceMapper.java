package org.rcsb.hadoop.io;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

/**
 * This class maps an encoded polymer chain from a Hadoop sequence file 
 * to a <PdbId ChainId, Polymer Sequence> pair.
 * 
 * @author Peter Rose
 */
public class HadoopToSequenceMapper implements PairFunction<Tuple2<Text,ArrayWritable>,String, String> {
	private static final long serialVersionUID = 1L;

	/**
	 * Maps an encoded polymer chain from a Hadoop sequence file to a <PdbId ChainId, Polymer Sequence> pair.
	 *
	 * @param  tuple a tuple of <PdbId.ChainID, Encoded Polymer Chain>
     * @return tuple with <PdbId.ChainID, Polymer Sequence>
     * 
	 * @see org.apache.spark.api.java.function.PairFunction#call(java.lang.Object)
	 */
	@Override
	public Tuple2<String, String> call(Tuple2<Text, ArrayWritable> tuple) throws Exception {
		Writable[] encodedPolymerChain = tuple._2.get();	
		SimplePolymerChain chain = new SimplePolymerChain(encodedPolymerChain);
		
		return new Tuple2<String, String>(tuple._1.toString(), chain.getSequence());
	}
}
