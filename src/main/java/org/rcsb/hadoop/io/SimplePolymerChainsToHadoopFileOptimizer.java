package org.rcsb.hadoop.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.vecmath.Point3d;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IOUtils;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.BZip2Codec;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.HTTPConnectionTools;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.rcsb.compress.AncientEgyptianDecomposition;
import org.rcsb.compress.CombinedTransform;
import org.rcsb.compress.DeltaTransform;
import org.rcsb.compress.FastWaveletTransform;
import org.rcsb.compress.IntegerDeltaZigzagVariableByte;
import org.rcsb.compress.IntegerToByteTransform;
import org.rcsb.compress.IntegerToByteTransformer;
import org.rcsb.compress.IntegerTransform;
import org.rcsb.compress.UnsignedDeltaTransform;
import org.rcsb.compress.dev.ByteQuadrupleToByteTransform;
import org.rcsb.compress.dev.ByteQuadrupleTransform;
import org.rcsb.compress.dev.DeltaHalfTransform;
import org.rcsb.compress.dev.DeltaToShortTransform;
import org.rcsb.compress.dev.PythagoreanQuadrupleTransform;
import org.rcsb.compress.dev.SphericalCoordinateTransform;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class creates a Hadoop sequence file for protein chains in the PDB. The Hadoop sequence file
 * uses a delta encoding of the PDB coordinates as well as BZIP2 block level compression.
 * 
 * 
 * @author  Peter Rose
 */
public class SimplePolymerChainsToHadoopFileOptimizer {
	private static final Logger logger = LoggerFactory.getLogger(SimplePolymerChainsToHadoopFileOptimizer.class);
	private static AtomCache cache = initializeCache();
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";
	

	public static void main(String[] args) throws Exception {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		String uri = args[0] 
				+ "_"
				+ timeStamp 
				+ ".seq";
		
		
		List<String> subset = new ArrayList<>(getAll());
//		List<String> pdbIds = subset;
		List<String> pdbIds = subset.subList(90000, 91000);

		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/peter/Data/PDB/");

		long start = System.nanoTime();
		// 82157
//		IntegerTransform t1 = new DeltaTransform();
		
		//
//		IntegerTransform t1 = new DeltaHalfTransform();
		
		// 82151
//		IntegerTransform t1 = new UnsignedDeltaTransform();
		
		
//		IntegerTransform t1 = new CombinedTransform(new DeltaTransform(), new PythagoreanQuadrupleTransform());
		
//		IntegerTransform t1 = new CombinedTransform(new DeltaTransform(), new ByteQuadrupleTransform());
		
		IntegerToByteTransform transform = new ByteQuadrupleToByteTransform(); // 1,321,934
//		IntegerToByteTransform transform = new DeltaToShortTransform();


		// 59157 ??? there are NaN in the output!!
//		IntegerTransform t0 = new CombinedTransform(new RefDeltaTransform(), new SphericalCoordinateTransform());
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 2")));

		// 59361
//		IntegerTransform t0 = new CombinedTransform(new RefDeltaTransform(), new SphericalCoordinateTransform());
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3")));
		
		// 71875
//		IntegerTransform t0 = new CombinedTransform(new DeltaTransform(), new SphericalCoordinateTransform());
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3")));

		// 72677
//		IntegerTransform t0 = new CombinedTransform(new NoOpTransform(), new SphericalCoordinateTransform());
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3")));
		
		// 74599
//		IntegerTransform t1 = new CombinedTransform(new DeltaTransform(), new SphericalCoordinateTransform()); 
		
		// 75486
//		IntegerTransform t0 = new CombinedTransform(new DeltaTransform(), new PythagoreanTransform()); 
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3"))); 

		// 76366
//		IntegerTransform t0 = new CombinedTransform(new DeltaTransform(), new PythagoreanTransform()); 
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 2"))); 
	
		// 77667
//		IntegerTransform t0 = new CombinedTransform(new DeltaTransform(), new ColToRowTransform());
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3")));

		// 77951
//		IntegerTransform t1 = new CombinedTransform(new DeltaTransform(), new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3"))); 

		// 90340
//		IntegerTransform t0 = new NoOpTransform(); 
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3"))); 
		
		// 93155
//		IntegerTransform t0 = new NoOpTransform(); 
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Daubechies 10"))); 
		
		//80575
//		IntegerTransform t0 = new RefDeltaTransform(); 
//		IntegerTransform t1 = new CombinedTransform(t0, new AncientEgyptianDecomposition(new FastWaveletTransform("Legendre 3"))); 

		// 94292
//		IntegerTransform t1 = new RefDeltaTransform();
		
		//		IntegerTransform t1 = new DeltaTransform();
		//		IntegerTransform t1 = new CombinedTransform(new NoOpTransform(), new AncientEgyptianDecomposition(new FastWaveletTransform("Daubechies 10")));
		
//		IntegerToByteTransform transform = new IntegerToByteTransformer(t1);
		
		toSequenceFile(uri, pdbIds, transform, true);
	
		long end = System.nanoTime();

		System.out.println("Time: " + (end - start)/1E9 + " sec.");
	}

	public static long toSequenceFile(String fileName, Collection<String> pdbIds, IntegerToByteTransform transform, boolean verbose)
			throws Exception {

		int failure = 0;
		int success = 0;
		int chains = 0;
		int[] metrics = new int[2];

		try (SequenceFile.Writer writer = SequenceFile.createWriter(new Configuration(),
				SequenceFile.Writer.file(new Path(fileName)),
				SequenceFile.Writer.keyClass(Text.class),
				SequenceFile.Writer.valueClass(SimplePolymerChain.class),
				SequenceFile.Writer.compression(SequenceFile.CompressionType.BLOCK, new BZip2Codec()));		
				) 
				{
			for (String pdbId: pdbIds) {
				logger.info(pdbId);

				Structure s = null;
				try {
					s = StructureIO.getStructure(pdbId);
					success++;
				} catch (Exception e) {
					// some files can't be read. Let's just skip those!
					logger.error(pdbId, e);
					failure++;
					continue;
				}

				if (s == null) {
					logger.error(pdbId + " structure object is null");
					continue;
				}

				if (s.getChains().size() == 0) {
					continue;
				}

				chains += append(writer, pdbId, s, transform, metrics);
			}
			IOUtils.closeStream(writer);
				}

		File f = new File(fileName);
		long length = f.length();
		
		if (verbose) {
			System.out.println("Total structures: " + pdbIds.size());
			System.out.println("Success: " + success);
			System.out.println("Failure: " + failure);
			System.out.println("File size: " + length);
			System.out.println("Chains: " + chains);
			System.out.println("Size: " + metrics[0]);
			System.out.println("Compression: "+  metrics[0]/(float) length);
			System.out.println("Time: " + metrics[1]/1E6);
		}
		
		return length;
	}

	private static int append(SequenceFile.Writer writer, String pdbId, Structure s, IntegerToByteTransform transform, int[] metrics)
			throws Exception {

		int chainCount = 0;

		for (Chain c: s.getChains()) {
			List<Group> groups = c.getAtomGroups(GroupType.AMINOACID);
			List<Group> nAcids = c.getAtomGroups(GroupType.NUCLEOTIDE);
//			List<Group> groups = c.getSeqResGroups(GroupType.AMINOACID);
//			List<Group> nAcids = c.getSeqResGroups(GroupType.NUCLEOTIDE);

			boolean aminoAcid = true;
			if (nAcids.size() > groups.size()) {
				groups = nAcids;
				aminoAcid = false;
			}
			int dna = 0;
			int rna = 0;
			int peptide = 0;
			int dPeptide = 0;
			int unknownResidues = 0;

			List<Point3d> coords = new ArrayList<>(groups.size());	

			StringBuilder sb = new StringBuilder();		

			int gaps = 0;

			for (int i = 0; i < groups.size(); i++) {
				Group g = groups.get(i);
				char code = g.getChemComp().getOne_letter_code().charAt(0);
				if (code == 'X') {
					unknownResidues++;
				}
				sb.append(code);

				PolymerType p = g.getChemComp().getPolymerType();
				if (p.equals(PolymerType.peptide)) {
					peptide++;
				} else if (p.equals(PolymerType.dpeptide)) {
					dPeptide++;
				} else if (p.equals(PolymerType.dna)) {
					dna++;
				} else if (p.equals(PolymerType.rna)) {
					rna++;
				}

				Atom atom = null;
				if (aminoAcid) {
					atom = ((AminoAcidImpl)g).getCA();
				} else {
					atom = ((NucleotideImpl)g).getP();
				}
				if (atom == null) {
					gaps++;
				} else {
					coords.add(new Point3d(atom.getCoords()));
				}
			}

			// ignore chains with less than 10 residues (with coordinates)
			//			System.out.println("size:  " + (groups.size()-gaps));
			if (groups.size()-gaps < 10) {
				continue;
			}

			if (unknownResidues > (groups.size()-gaps)/2) {
				logger.info(pdbId + ": Polymer with many unknown residues ignored: " + pdbId + c.getChainID());
				continue;
			}
			// ignore any mixed polymer types
			if (dPeptide > 0) {
				logger.info(pdbId  + c.getChainID() + ": d-peptide ignored");
				continue;
			}
			if (dna > 0 && rna > 0) {
				logger.info(pdbId + c.getChainID() + ": DNA/RNA hybrid ignored");
				continue;
			}

			// determine polymer type
			SimplePolymerType polymerType = null;
			if (peptide > 0) {
				polymerType = SimplePolymerType.PROTEIN;
			} else if (dna > 0) {
				polymerType = SimplePolymerType.DNA;
			} else if (rna > 0) {
				polymerType = SimplePolymerType.RNA;
			} else {
				continue;
			}

			Point3d[] coordinates = coords.toArray(new Point3d[coords.size()]);
			chainCount++;

			Text key = new Text(pdbId + "." + c.getChainID());
			System.out.println(key);
			
			
			long start = System.nanoTime();
			SimplePolymerChain value = new SimplePolymerChain(transform);
//			SimplePolymerChain value = new SimplePolymerChain();
			value.setPolymerType(polymerType.ordinal());
			value.setCoordinates(coordinates);
			value.setSequence(sb.toString());
		    metrics[1] += System.nanoTime() - start;
			
			writer.append(key, value);
			
			metrics[0] += value.size();
		}
		return chainCount;
	}

	/**
	 * Returns the current list of all PDB IDs.
	 * @return PdbChainKey set of all PDB IDs.
	 */
	public static SortedSet<String> getAll() {
		SortedSet<String> representatives = new TreeSet<String>();

		try {

			URL u = new URL(allUrl);

			InputStream stream = HTTPConnectionTools.getInputStream(u, 60000);

			if (stream != null) {
				BufferedReader reader = new BufferedReader(
						new InputStreamReader(stream));

				String line = null;

				while ((line = reader.readLine()) != null) {
					int index = line.lastIndexOf("structureId=");
					if (index > 0) {
						representatives.add(line.substring(index + 13, index + 17));
					}
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		return representatives;
	}

	private static AtomCache initializeCache() {
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
//		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
//		params.setParseCAOnly(true); // can't use CA only since we need to read DNA/RNA
//		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(false);
		cache.setFileParsingParams(params);
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
		return cache;
	}
}
