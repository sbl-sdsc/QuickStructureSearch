package org.rcsb.hadoop.io;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

/**
 * This class maps an encoded polymer chain from a Haddop sequence file to a <PdbId.ChainId, Polymer Coordinates> pair.
 *
 * @author Peter Rose
 *
 */
public class HadoopToCoordinateMapper implements PairFunction<Tuple2<Text,ArrayWritable>,String, Point3d[]> {
	private static final long serialVersionUID = 1L;

	/**
	 * Maps an encoded polymer chain from a Haddop sequence file to a <PdbId ChainId, Point3d[]> pair.
	 * Note, the Point3d[] array contains null entries at gaps in the protein chain.
	 * 
	 * @param  tuple a tuple of <PdbId.ChainID, Encoded Polymer Chain>
     * @return tuple with <PdbId.ChainID, Polymer Coordinates>
     * 
	 * @see org.apache.spark.api.java.function.PairFunction#call(java.lang.Object)
	 */
	@Override
	public Tuple2<String, Point3d[]> call(Tuple2<Text, ArrayWritable> tuple) throws Exception {
		Writable[] encodedPolymerChain = tuple._2.get();
        SimplePolymerChain chain = new SimplePolymerChain(encodedPolymerChain);
        
		return new Tuple2<String,Point3d[]>(tuple._1.toString(), chain.getCoordinates());
	}
}
