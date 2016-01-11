package org.rcsb.hadoop.io;

import org.apache.hadoop.io.Text;
import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;
/**
 * This class maps a Tuple2<Text,SimplePolymerChain> read from a Hadoop sequence files
 * to a Tuple2<String, SimplePolymerChain> that can be cached in Spark. Note, the original tuple read
 * from a Hadoop sequence file cannot be cached, since the Text and SimplePolymerChain objects are
 * reused by the Hadoop sequence file reader. This class converts the Text object to String and makes
 * a copy of the SimplePolymerChain object.
 *
 * @author Peter Rose
 *
 */
public class HadoopToSimplePolymerChainMapper implements PairFunction<Tuple2<Text,SimplePolymerChain>,String, SimplePolymerChain> {
	private static final long serialVersionUID = 1L;

	/**
	 * Maps a Tuple2<Text,SimplePolymerChain> to a Tuple2<Text,SimplePolymerChain>
	 * 
	 * @param  tuple a tuple of <Text, SimplePolymerChain>
     * @return tuple q tuple of <String, SimplePolymerChain>
     * 
	 * @see org.apache.spark.api.java.function.PairFunction#call(java.lang.Object)
	 */
	@Override
	public Tuple2<String, SimplePolymerChain> call(Tuple2<Text, SimplePolymerChain> t) throws Exception {
		return new Tuple2<String, SimplePolymerChain>(t._1.toString(), new SimplePolymerChain(t._2));
	}
}
