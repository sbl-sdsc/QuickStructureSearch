package org.rcsb.structuralSimilarity;

import java.io.IOException;

import javax.vecmath.Point3d;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.io.compress.BZip2Codec;

public class PdbChainToSequenceFileWriter {
	private static final int SCALE = 1000; // Factor to convert PDB coordinates to integer values. Do not change this value!
	private SequenceFile.Writer writer;
	private Text key = new Text();
	private ArrayWritable value = new IntArrayWritable();	

	public PdbChainToSequenceFileWriter(String uri) throws IOException {
		writer = SequenceFile.createWriter(new Configuration(),
				SequenceFile.Writer.file(new Path(uri)),
				SequenceFile.Writer.keyClass(key.getClass()),
				SequenceFile.Writer.valueClass(value.getClass()),
				SequenceFile.Writer.compression(SequenceFile.CompressionType.BLOCK, new BZip2Codec()));	
	}

	public void write(String chainId, Point3d[] points) throws IOException {
		int gaps = 0;
		for (Point3d p: points) {
			if (p == null) {
				gaps++;
			}
		}

		Writable[] writable = new Writable[(points.length+1)*3 - 2*gaps+1];
		int n = 0;

		// record number of points
		writable[n++] = new IntWritable(points.length);

		int x = 0;
		int y = 0;
		int z = 0;

		// delta encode coordinate values. Gaps in the protein chains
		// are encoded as the maximum integer values.
		for (int i = 0, dx = 0, dy = 0, dz = 0; i < points.length; i++) {
			Point3d p = points[i];
			if (p != null) {
				// delta encode coordinate values as integers
				x = (int)Math.round(p.x*SCALE);
				y = (int)Math.round(p.y*SCALE);
				z = (int)Math.round(p.z*SCALE);
				writable[n++] = new IntWritable(x-dx);
				writable[n++] = new IntWritable(y-dy);
				writable[n++] = new IntWritable(z-dz);
				dx = x;
				dy = y;
				dz = z;
			} else {
				// encode a gap in the protein chain
				writable[n++] = new IntWritable(Integer.MAX_VALUE);
			}
		}

		// record last x,y,z values for validation
		writable[n++] = new IntWritable(x);
		writable[n++] = new IntWritable(y);
		writable[n++] = new IntWritable(z);

		key.set(chainId);
		value.set(writable);
		writer.append(key, value);
	}

	public void close() throws IOException {
		writer.close();
	}
}
