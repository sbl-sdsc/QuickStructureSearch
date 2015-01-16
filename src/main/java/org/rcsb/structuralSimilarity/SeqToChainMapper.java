package org.rcsb.structuralSimilarity;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

public class SeqToChainMapper implements PairFunction<Tuple2<Text,ArrayWritable>,String, Point3d[]> {
	private static final long serialVersionUID = 1L;
	private static final double SCALE = 0.001;

	public SeqToChainMapper() {
	}

	@Override
	public Tuple2<String, Point3d[]> call(Tuple2<Text, ArrayWritable> tuple) throws Exception {
		Writable[] w = tuple._2.get();
		int len = ((IntWritable)w[0]).get();
		Point3d[] points = new Point3d[len];
		
		int j = 1;
		int x = 0;
		int y = 0;
		int z = 0;

		for (int i = 0; i < points.length; i++) {
			int v = ((IntWritable)w[j++]).get();
			if (v == Integer.MAX_VALUE) {
				points[i] = null; // a gap in the coordinates is represented by a null value
			} else {
				x += v;
				y += ((IntWritable)w[j++]).get();
				z += ((IntWritable)w[j++]).get();
				points[i] = new Point3d(x*SCALE, y*SCALE, z*SCALE);
			}
		}
		
		// compare the last x, y, z values with the expected values
		if (x != ((IntWritable)w[j++]).get()) {
			throw new Exception("ERROR: Input file is corrupted");
		}
		if (y != ((IntWritable)w[j++]).get()) {
			throw new Exception("ERROR: Input file is corrupted");
		}
		if (z != ((IntWritable)w[j++]).get()) {
			throw new Exception("ERROR: Input file is corrupted");
		}
		
		return new Tuple2<String,Point3d[]>(tuple._1.toString(), points);

	}
}
