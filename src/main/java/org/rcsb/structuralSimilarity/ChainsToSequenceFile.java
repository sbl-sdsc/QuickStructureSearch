package org.rcsb.structuralSimilarity;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.IOUtils;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.io.compress.BZip2Codec;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.rcsb.utils.BlastClustReader;
/**
 * This class creates a Hadoop sequence file for protein chains. By default, it uses the first
 * representative chains from each 40% sequence identity cluster. The Hadoop sequence file
 * uses a delta encoding of the PDB coordinates as well as BZIP2 block level compression.
 * 
 * @author  Peter Rose
 */
public class ChainsToSequenceFile {
	private static AtomCache cache = initializeCache();
	private static final int SCALE = 1000; // Factor to convert PDB coordinates to integer values. Do not change this value!

	private static final int PERCENT_SEQUENCE_IDENTITY = 40; // use 40% sequence identity clusters

	public static void main(String[] args) throws IOException {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		
		String uri = "/Users/peter/Data/PDB_CHAINS/protein_chains_" 
				+ PERCENT_SEQUENCE_IDENTITY 
				+ "_"
				+ timeStamp 
				+ ".seq";

		BlastClustReader reader = new BlastClustReader(PERCENT_SEQUENCE_IDENTITY);
		List<List<String>> clusters = reader.getPdbChainIdClusters();
		clusters = clusters.subList(0, 250);
	
		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/peter/Data/PDB/");
		
		long start = System.nanoTime();
		
		Text key = new Text();
		ArrayWritable value = new IntArrayWritable();	
		SequenceFile.Writer writer = null;

		int failure = 0;
		int success = 0;
		int chains = 0;
		
		try {
			writer = SequenceFile.createWriter(new Configuration(),
					SequenceFile.Writer.file(new Path(uri)),
					SequenceFile.Writer.keyClass(key.getClass()),
					SequenceFile.Writer.valueClass(value.getClass()),
			        SequenceFile.Writer.compression(SequenceFile.CompressionType.BLOCK, new BZip2Codec()));	

			for (List<String> cluster: clusters) {
				System.out.println(cluster.get(0));
				Structure s = null;
				try {
					s = StructureIO.getStructure(cluster.get(0));
					success++;
				} catch (Exception e) {
					// some files can't be read. Let's just skip those!
					e.printStackTrace();
					failure++;
					continue;
				}
				
				if (s == null) {
					System.err.println("structure null: " + cluster.get(0));
					continue;
				}

				if (s.getChains().size() == 0) {
					continue;
				}
				
				int gaps = 0;
				Chain c = s.getChain(0);
				List<Group> groups = c.getSeqResGroups(GroupType.AMINOACID);
				// does this list include alternate locations?
				List<Atom> cas = new ArrayList<>(groups.size());
				for (Group g: groups) {
					if (g instanceof AminoAcidImpl) {
						Atom ca = ((AminoAcidImpl)g).getCA();
						cas.add(ca);
						if (ca == null) {
							gaps++;
						}
					}
				}
				// ignore chains with less than 24 residues (with coordinates)
				if (cas.size()-gaps < 24) {
					continue;
				}
				chains++;
				if (chains % 1000 == 0) {
					System.out.println("Chains saved: " + chains);
				}

				key.set(cluster.get(0));
				value.set(toWritable(cas.toArray(new Atom[0]), gaps));
				writer.append(key, value);
			}
		} finally {
			IOUtils.closeStream(writer);
		}
		
		System.out.println("Total structures: " + clusters.size());
		System.out.println("Success: " + success);
		System.out.println("Chains: " + chains);
		System.out.println("Failure: " + failure);
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	private static Writable[] toWritable(Atom[] ca, int gaps) {
		Writable[] writable = new Writable[(ca.length+1)*3 - 2*gaps + 1];
		int n = 0;
		
		// record number of points
		writable[n++] = new IntWritable(ca.length);
		
		int x = 0;
		int y = 0;
		int z = 0;

		// delta encode coordinate values. Gaps in the protein chains
		// are encoded as the maximum integer values.
		for (int i = 0, dx = 0, dy = 0, dz = 0; i < ca.length; i++) {
			Atom a = ca[i];
			if (a != null) {
				// delta encode coordinate values as integers
				x = (int)Math.round(a.getX()*SCALE);
				y = (int)Math.round(a.getY()*SCALE);
				z = (int)Math.round(a.getZ()*SCALE);
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
		
		return writable;
	}

	private static AtomCache initializeCache() {
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
//		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
//		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(false);
		cache.setFileParsingParams(params);
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
//		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		return cache;
	}
}
