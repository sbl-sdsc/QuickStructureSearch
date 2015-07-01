package org.rcsb.hadoop.io;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;
/**
 * This class maps an encoded polymer chain Tuple from a Hadoop sequence file 
 * to a Tuple of <PdbId ChainId, CompactPolymerChain>.
 *
 * @author Peter Rose
 *
 */
public class HadoopToSimpleChainMapper implements PairFunction<Tuple2<Text,ArrayWritable>,String, SimplePolymerChain> {
	private static final long serialVersionUID = 1L;

	/**
	 * Maps an encoded polymer chain <Text, ArrayWritable> pair to a <PdbId.chainId, SimplePolymerChain> pair.
	 * 
	 * @param  tuple a tuple of <PdbId.ChainID, Encoded Polymer Chain>
     * @return tuple q tuple of <PdbId.ChainID, SimplePolymerChain>
     * 
	 * @see org.apache.spark.api.java.function.PairFunction#call(java.lang.Object)
	 */
	@Override
	public Tuple2<String, SimplePolymerChain> call(Tuple2<Text, ArrayWritable> tuple) throws Exception {
		Writable[] encodedPolymerChain = tuple._2.get();
		SimplePolymerChain chain = new SimplePolymerChain(encodedPolymerChain);
		return new Tuple2<String, SimplePolymerChain>(tuple._1.toString(), chain);
	}
}
