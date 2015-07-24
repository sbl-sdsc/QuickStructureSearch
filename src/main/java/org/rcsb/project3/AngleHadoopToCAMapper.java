package org.rcsb.project3;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFunction;
import org.rcsb.hadoop.io.SimplePolymerChain;

import scala.Tuple2;

/**
 * This class maps an encoded polymer chain from a Haddop sequence file to a <PdbId.ChainId, Polymer Coordinates> pair.
 *
 * @author Peter Rose
 *
 */
public class AngleHadoopToCAMapper implements PairFunction<Tuple2<String, SimplePolymerChain>,String, Point3d[]> {
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
	public Tuple2<String, Point3d[]> call(Tuple2<String, SimplePolymerChain> tuple) throws Exception {
		Point3d[] atoms = tuple._2.getCoordinates();
		Point3d[] CA = new Point3d[atoms.length/3];
		for (int i = 0, j = 0; i < atoms.length; i+=3, j++) {
			CA[j] = atoms[i];
		}
		return new Tuple2<String,Point3d[]>(tuple._1.toString(), CA);
	}
}
