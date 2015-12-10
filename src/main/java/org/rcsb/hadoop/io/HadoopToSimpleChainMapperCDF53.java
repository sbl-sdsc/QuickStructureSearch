package org.rcsb.hadoop.io;

import java.io.Serializable;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.spark.api.java.function.PairFunction;
import org.rcsb.compress.IntegerTransform;

import scala.Tuple2;
/**
 * This class maps an encoded polymer chain Tuple from a Hadoop sequence file 
 * to a Tuple of <PdbId ChainId, CompactPolymerChain>.
 *
 * @author Peter Rose
 *
 */
public class HadoopToSimpleChainMapperCDF53 implements PairFunction<Tuple2<Text,ArrayWritable>,String, SimplePolymerChainCDF53>, Serializable {
	private static final long serialVersionUID = 1L;
	private IntegerTransform transform;

	public HadoopToSimpleChainMapperCDF53(IntegerTransform transform) {
		this.transform = transform;
	}
	
	/**
	 * Maps an encoded polymer chain <Text, ArrayWritable> pair to a <PdbId.chainId, SimplePolymerChain> pair.
	 * 
	 * @param  tuple a tuple of <PdbId.ChainID, Encoded Polymer Chain>
     * @return tuple q tuple of <PdbId.ChainID, SimplePolymerChain>
     * 
	 * @see org.apache.spark.api.java.function.PairFunction#call(java.lang.Object)
	 */
	@Override
	public Tuple2<String, SimplePolymerChainCDF53> call(Tuple2<Text, ArrayWritable> tuple) throws Exception {
		Writable[] encodedPolymerChain = tuple._2.get();
		SimplePolymerChainCDF53 chain = new SimplePolymerChainCDF53(encodedPolymerChain, transform);
		return new Tuple2<String, SimplePolymerChainCDF53>(tuple._1.toString(), chain);
	}
}
