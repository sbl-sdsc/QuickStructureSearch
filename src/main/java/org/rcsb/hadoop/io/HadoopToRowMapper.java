package org.rcsb.hadoop.io;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.RowFactory;

import scala.Tuple2;
/**
 * This class maps an encoded polymer chain Tuple from a Hadoop sequence file 
 * to a Row object
 *
 * @author Peter Rose
 *
 */
public class HadoopToRowMapper implements Function<Tuple2<Text,ArrayWritable>, Row> {
	private static final long serialVersionUID = 1L;

	/**
	 * Maps an encoded polymer chain <Text, ArrayWritable> pair to a Row.
	 * 
	 * @param  tuple a tuple of <PdbId.ChainID, Encoded Polymer Chain>
     * @return row a row in a DataFrame
     * 
	 */
	@Override
	public Row call(Tuple2<Text, ArrayWritable> tuple) throws Exception {
		String key = tuple._1.toString();
		String hashCol = key.substring(1,2);
		Writable[] encodedPolymerChain = tuple._2.get();	
		int[] values = new int[encodedPolymerChain.length];
		
		for (int i = 0; i < encodedPolymerChain.length; i++) {
			values[i] = ((IntWritable)encodedPolymerChain[i]).get();
		}
		return RowFactory.create(hashCol,key, values);
	}
}
