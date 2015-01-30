package org.rcsb.structuralSimilarity;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;

import scala.Tuple2;

/**
 * Demo program to show how to smooth protein chains and save them to pdb files for
 * visualization.
 * To Visualize a .pdb file, launch the Protein Workshop viewer using this URL
 * http://www.rcsb.org/pdb/explore/viewerLaunch.do?viewerType=PW&structureId=1STP&unit=bio&unit_id=1&useTrace=n
 * Then use the File->Open File menu to read a pdb file
 * @author  Peter Rose
 */
public class SmoothTest { 
	private static int NUM_THREADS = 8;

	public static void main(String[] args ) throws FileNotFoundException
	{
		String sequenceFileName = "src/test/resources/protein_chains_40_20150114_141156.seq";

		SmoothTest aaa = new SmoothTest();
		aaa.run(sequenceFileName, args[0]);
	}

	private void run(String path, String filename) throws FileNotFoundException {
		// setup spark
		SparkConf conf = new SparkConf()
		.setMaster("local[" + NUM_THREADS + "]")
		.setAppName("1" + this.getClass().getSimpleName())
		.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);

		// Step 1. calculate <pdbId.chainId, feature vector> pairs
		List<Tuple2<String,Point3d[]>> list  = sc
				.sequenceFile(path, Text.class, ArrayWritable.class, NUM_THREADS)  // read protein chains
				.sample(false, 0.001, 123456) // use only a random fraction, i.e., 40%
				.mapToPair(new SeqToChainMapper()) // convert input to <pdbId.chainId, CA coordinate> pairs
				.filter(new GapFilter(3, 5)) // keep protein chains with gap size <= 3 and <= 5 gaps
				.filter(new LengthFilter(50,1000)) // keep protein chains with at least 75 residues
		//		.mapToPair(new ChainSmootherMapper(new RogenSmoother(2))); // apply smoothing here ..
				//		.mapToPair(new ChainSmootherMapper(new SavitzkyGolay7PointSmoother(2))); // apply smoothing here ..
				.collect(); // return results to master node

		writePointsToPdb(list, filename);

		sc.stop();
		sc.close();
	}

	/**
	 * Writes each tuple of <PdbID.ChainID, Point3[]> as a separate .pdb file
	 * @param smoothed
	 * @param filename
	 * @throws FileNotFoundException
	 */
	private static void writePointsToPdb(List<Tuple2<String,Point3d[]>> smoothed, String fileName) throws FileNotFoundException {	
		for (Tuple2<String,Point3d[]> tuple: smoothed) {
			PrintWriter writer = new PrintWriter(fileName + tuple._1 + ".pdb");
			
			// create a protein chain object
			Chain c = new ChainImpl();
			c.setChainID(tuple._1.substring(5)); // this substring represents the chain id
			
			Point3d[] points = tuple._2;
			for (int i = 0; i < points.length; i++) {
				if (points[i] != null) {
					// create a C-alpha atom
					Atom atom = new AtomImpl();
					atom.setName("CA");
					atom.setPDBserial(i);
					atom.setAltLoc(' ');
					atom.setX(points[i].x);
					atom.setY(points[i].y);
					atom.setZ(points[i].z);
					
					// create an amino acid group (residue) and add it to the chain
					Group g = new AminoAcidImpl();
					g.setPDBName("GLY"); // for now, all amino acids set to Glycine
					g.setResidueNumber(tuple._1.substring(5,6),  i, ' ');
					g.addAtom(atom);
					
					// add group to chain
					c.addGroup(g);
					
                    // write atom to PDB file
					writer.println(atom.toPDB());
				}
			}
			writer.close();
		}
	}
}
