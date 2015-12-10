package org.rcsb.hadoop.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.vecmath.Point3d;

import me.lemire.integercompression.BinaryPacking;
import me.lemire.integercompression.Composition;
import me.lemire.integercompression.DeltaZigzagBinaryPacking;
import me.lemire.integercompression.DeltaZigzagVariableByte;
import me.lemire.integercompression.FastPFOR;
import me.lemire.integercompression.FastPFOR128;
import me.lemire.integercompression.IntCompressor;
import me.lemire.integercompression.IntegerCODEC;
import me.lemire.integercompression.Kamikaze;
import me.lemire.integercompression.NewPFD;
import me.lemire.integercompression.NewPFDS16;
import me.lemire.integercompression.OptPFDS16;
import me.lemire.integercompression.Simple16;
import me.lemire.integercompression.SkippableComposition;
import me.lemire.integercompression.SkippableIntegerCODEC;
import me.lemire.integercompression.VariableByte;
import me.lemire.integercompression.differential.IntegratedIntCompressor;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.IOUtils;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.BZip2Codec;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.io.compress.Lz4Codec;
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
import org.rcsb.compress.Delta4Transform;
import org.rcsb.compress.DeltaReverseTransform;
import org.rcsb.compress.DeltaTransform;
import org.rcsb.compress.HarrWavelet;
import org.rcsb.compress.HorizontalDeltaTransform;
import org.rcsb.compress.IntCompressorTransform;
import org.rcsb.compress.IntegerTransform;
import org.rcsb.compress.LeGallWavelet;
import org.rcsb.compress.NullOpTransform;
import org.rcsb.compress.PFORTransform;
import org.rcsb.compress.PythagoreanTransform;
import org.rcsb.compress.UnsignedDeltaTransform;
import org.rcsb.compress.UnsignedPythagoreanTransform;
import org.rcsb.compress.UnsignedRefDeltaTransform;
import org.rcsb.compress.UnsignedTransform;
import org.rcsb.compress.XorTransform;
import org.rcsb.structuralSimilarity.IntArrayWritable;

/**
 * This class creates a Hadoop sequence file for protein chains in the PDB. The Hadoop sequence file
 * uses a delta encoding of the PDB coordinates as well as BZIP2 block level compression.
 * 
 * Example run:
 * Jan. 20 PDB release
 * Total structures: 105906
 * Success: 105899
 * Failure: 7
 * Chains: 280565
 * Time: 50605.43170881 sec.
 * Size: 341 MB
 * 
 * @author  Peter Rose
 */
public class NewPdbChainsToHadoopFile {
	private static AtomCache cache = initializeCache();
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";

	public static void main(String[] args) throws IOException {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

//		String uri = args[0] 
//				+ "_"
//				+ timeStamp 
//				+ ".seq";
		
		String uri = args[0] 
				+ "_"
				+ "latest" 
				+ ".seq";

		//		List<String> dna = Arrays.asList("4RNK","4RO4","4RO7");
//		List<String> protein = Arrays.asList("1STP");
//		List<String> pdbIds = protein;
//		Set<String> pdbIds = getAll();

		// create a subset
		List<String> subset = new ArrayList<>(getAll());
//		List<String> pdbIds = subset.subList(0, 100);
		List<String> pdbIds = subset.subList(0, 100);


		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/peter/Data/PDB/");

		long start = System.nanoTime();
		toSequenceFile(uri, pdbIds, true);
		long end = System.nanoTime();

		System.out.println("Time: " + (end - start)/1E9 + " sec.");
	}

	public static void toSequenceFile(String fileName, Collection<String> pdbIds, boolean verbose)
			throws IOException {

		int failure = 0;
		int success = 0;
		int chains = 0;
		int[] metrics = new int[2];
		long start = 0;
		long end = 0;
		long time = 0;

		try (SequenceFile.Writer writer = SequenceFile.createWriter(new Configuration(),
				SequenceFile.Writer.file(new Path(fileName)),
				SequenceFile.Writer.keyClass(Text.class),
				SequenceFile.Writer.valueClass(IntArrayWritable.class),
				SequenceFile.Writer.compression(SequenceFile.CompressionType.BLOCK, new BZip2Codec()));		
				) 
				{
			for (String pdbId: pdbIds) {
				if (verbose) {
					System.out.println(pdbId);
				}

				Structure s = null;
				try {
					s = StructureIO.getStructure(pdbId);
					success++;
				} catch (Exception e) {
					// some files can't be read. Let's just skip those!
					e.printStackTrace();
					failure++;
					continue;
				}

				if (s == null) {
					System.err.println("structure null: " + pdbId);
					continue;
				}

				if (s.getChains().size() == 0) {
					continue;
				}

				start = System.nanoTime();
				chains += append(writer, pdbId, s, metrics);
				time += System.nanoTime() - start;
			}
			IOUtils.closeStream(writer);
				}

		if (verbose) {
			System.out.println("Total structures: " + pdbIds.size());
			System.out.println("Success: " + success);
			System.out.println("Failure: " + failure);
			System.out.println("Chains: " + chains);
			System.out.println("Size: " + metrics[0]);
			System.out.println("Time: " + metrics[1]/1E6);
		}
	}

	private static int append(SequenceFile.Writer writer, String pdbId, Structure s, int[] metrics)
			throws IOException {

//       IntegerTransform transform = new IntCompressorTransform(new SkippableComposition(new BinaryPacking(), new VariableByte())); // 109K size: 35045
//		IntegerTransform transform = new IntCompressorTransform(new SkippableComposition(new FastPFOR(), new VariableByte())); // 110K size: 33520 Time: 62.954806 (forward Time: 30.78215)
		IntegerTransform transform = new PFORTransform(new Composition(new DeltaZigzagVariableByte(), new VariableByte())); // 80K size: 26548 Time: 44.758336 Time: 54.547389 (forward: Time: 19.477541)
//		IntegerTransform transform = new PFORTransform(new Composition(new FastPFOR(), new DeltaZigzagVariableByte())); // 100K size: 30175
//		IntegerTransform transform = new IntCompressorTransform(new SkippableComposition(new Simple16(), new VariableByte())); // fail
//		IntegerTransform transform = new IntCompressorTransform(new SkippableComposition(new NewPFD(), new VariableByte())); // fail
//		IntegerTransform transform = new IntCompressorTransform(new Simple16()); // fail
//		IntegerTransform transform = new IntCompressorTransform(new NewPFD()); // fail
//		IntegerTransform transform = new IntCompressorTransform(new VariableByte()); // fail
//		IntegerTransform transform = new IntCompressorTransform(new BinaryPacking()); // fail
//		IntegerTransform transform = new IntCompressorTransform(new FastPFOR()); // fail
    //    IntegerTransform transform = new IntCompressorTransform(new Kamikaze()); // fail
//		IntegerTransform transform = new DeltaTransform(); // 80K, size: 59613
//		IntegerTransform transform = new HorizontalDeltaTransform(); //
//		IntegerTransform transform = new NullOpTransform(); // 100K size: 59613
//		IntegerTransform transform = new UnsignedDeltaTransform(); // 80K
//		IntegerTransform transform = new UnsignedRefDeltaTransform(); // 80K
//		IntegerTransform transform = new XorTransform(); // 88K
//		IntegerTransform transform = new CombinedTransform(new DeltaTransform(), new HorizontalDeltaTransform()); // 84K
//		IntegerTransform transform1 = new CombinedTransform(new DeltaTransform(), new HorizontalDeltaTransform());
//		IntegerTransform transform = new CombinedTransform(transform1, new PFORTransform(new DeltaZigzagVariableByte())); // 88K
//		IntegerTransform transform = new PFORTransform(new DeltaZigzagVariableByte()); // 80K size: 26548
//		IntegerTransform transform = new PFORTransform(new Composition(new NewPFDS16(),new VariableByte())); // 113K Size: 35016
//		IntegerTransform transform = new PFORTransform(new Composition(new Kamikaze(),new VariableByte())); // fail
//		IntegerTransform transform = new AncientEgyptianDecomposition(new LeGallWavelet()); // 88K
//		IntegerTransform transform = new AncientEgyptianDecomposition(new HarrWavelet()); // 105K
//		IntegerTransform transform = new CombinedTransform(new DeltaTransform(), new AncientEgyptianDecomposition(new HarrWavelet())); // 95K
//		IntegerTransform transform1 = new AncientEgyptianDecomposition(new HarrWavelet()); 
//		IntegerTransform transform = new CombinedTransform(transform1, new PFORTransform(new DeltaZigzagVariableByte())); // 103 Size: 29882
//		IntegerTransform transform = new CombinedTransform(new AncientEgyptianDecomposition(new LeGallWavelet()), new PFORTransform(new DeltaZigzagVariableByte())); // 91K Size: 27645 Time: 57.734268
//		IntegerTransform transform = new AncientEgyptianDecomposition(new PFORTransform(new DeltaZigzagVariableByte())); // index out of bounds
//		IntegerTransform transform = new AncientEgyptianDecomposition(new DeltaTransform()); // 81K
//		IntegerTransform transform = new AncientEgyptianDecomposition(new PFORTransform(new FastPFOR())); // fail, arrays need to be the same length
//		IntegerTransform transform = new CombinedTransform(new DeltaTransform(), new PythagoreanTransform()); // reverse not implemented // 78K
//		IntegerTransform transform1 = new CombinedTransform(new DeltaTransform(), new PythagoreanTransform());
//		IntegerTransform transform = new CombinedTransform(transform1, new PFORTransform(new DeltaZigzagVariableByte())); // size: 27274
//		IntegerTransform transform = new CombinedTransform(new UnsignedDeltaTransform(), new PFORTransform(new Composition(new FastPFOR(),new VariableByte()))); // 86K Size: 24441
//		IntegerTransform transform = new CombinedTransform(new DeltaTransform(), new PFORTransform(new Composition(new FastPFOR(),new VariableByte())));
//		IntegerTransform transform = new UnsignedTransform();
//		IntegerTransform transform = new PFORTransform(new DeltaZigzagVariableByte());
//		IntegerTransform transform = new PFORTransform(new DeltaZigzagBinaryPacking()); // no match?
//		IntegerTransform transform = new CombinedTransform(new UnsignedDeltaTransform(),new AncientEgyptianDecomposition(new LeGallWavelet()));
//		IntegerTransform transform = new CombinedTransform(new UnsignedDeltaTransform(), new PFORTransform(new FastPFOR()));
//		IntegerTransform transform = new CombinedTransform(new UnsignedDeltaTransform(), new PFORTransform(new FastPFOR()));
//		IntegerTransform transform = new CombinedTransform(new UnsignedTransform(), new PFORTransform(new FastPFOR()));
//		IntegerTransform transform1 = new CombinedTransform(new DeltaTransform(), new PythagoreanTransform());

//		IntegerTransform transform1 = new CombinedTransform(new DeltaTransform(), new UnsignedPythagoreanTransform());
//		IntegerTransform transform = new CombinedTransform(transform1, new PFORTransform(new FastPFOR()));
//		IntegerTransform transform = new Delta4Transform();
//		IntegerTransform transform = new DeltaReverseTransform();
//		IntegerTransform transform = new CombinedTransform(new AncientEgyptianDecomposition(new LeGallWavelet()), new PFORTransform(new Simple16()));
//		IntegerTransform transform = new PFORTransform(new VariableByte());
//		IntegerTransform transform = new UnsignedDeltaTransform();
//		IntegerTransform transform = new PFORTransform(new Composition(new FastPFOR(),new VariableByte()));
//		IntegerTransform transform1 = new PFORTransform(new VariableByte());
//		IntegerTransform transform = new CombinedTransform(new UnsignedDeltaTransform(), transform1);
//		IntegerTransform transform = new CombinedTransform(new DeltaTransform(), new PFORTransform(new Composition(new NewPFDS16(),new VariableByte()))); // 90K Size: 58506
//		IntegerTransform transform = new AncientEgyptianDecomposition(new LeGallWavelet());
//		IntegerTransform transform1 = new AncientEgyptianDecomposition(new LeGallWavelet());
//		IntegerTransform transform = new CombinedTransform(transform1, new PFORTransform(new DeltaZigzagVariableByte())); // 91K Size: 27645
//		IntegerTransform transform = new CombinedTransform(new DeltaTransform(), new AncientEgyptianDecomposition(new LeGallWavelet()));
		
		int chainCount = 0;

		for (Chain c: s.getChains()) {
			List<Group> groups = c.getSeqResGroups(GroupType.AMINOACID);
//			System.out.println("seq len: " + c.getSeqResSequence().length());
//			System.out.println("seqresgroups.size(): " + c.getSeqResGroups(GroupType.AMINOACID).size());
			List<Group> nAcids = c.getSeqResGroups(GroupType.NUCLEOTIDE);

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

			Point3d[] coords = new Point3d[groups.size()];	
	//		Integer[] sequence = new Integer[groups.size()];
			StringBuilder sb = new StringBuilder();		

			int gaps = 0;

			for (int i = 0; i < groups.size(); i++) {
				Group g = groups.get(i);
				char code = g.getChemComp().getOne_letter_code().charAt(0);
				if (code == 'X') {
					unknownResidues++;
				}
				sb.append(code);
	//			sequence[i] = (int)code;

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

//				for (Atom a: g.getAtoms()) {
//					System.out.println(a.toPDB());
//				}
				Atom atom = null;
				if (aminoAcid) {
					atom = ((AminoAcidImpl)g).getCA();
				} else {
					atom = ((NucleotideImpl)g).getP();
				}
				if (atom == null) {
					gaps++;
				} else {
					coords[i] = new Point3d(atom.getCoords());
				}
			}

			// ignore chains with less than 10 residues (with coordinates)
			//			System.out.println("size:  " + (groups.size()-gaps));
			if (groups.size()-gaps < 10) {
				continue;
			}

			if (unknownResidues > (groups.size()-gaps)/2) {
				System.err.println("Polymer with many unknown residues ignored: " + pdbId + c.getChainID());
				continue;
			}
			// ignore any mixed polymer types
			if (dPeptide > 0) {
				System.err.println("d-peptide ignored: " + pdbId + c.getChainID());
				continue;
			}
			if (dna > 0 && rna > 0) {
				System.err.println("DNA/RNA hybrid ignored: " + pdbId + c.getChainID());
				continue;
			}

			// determine polymer type
			SimplePolymerType polymerType = null;
			if (peptide > 0) {
				polymerType = SimplePolymerType.PROTEIN;
				//				System.out.println("PROTEIN: " + dna);
			} else if (dna > 0) {
				polymerType = SimplePolymerType.DNA;
				//				System.out.println("DNA: " + dna);
			} else if (rna > 0) {
				polymerType = SimplePolymerType.RNA;
				//				System.out.println("RNA: " + rna);
			} else {
				continue;
			}

			chainCount++;

			Text key1 = new Text(pdbId + "." + c.getChainID());
			ArrayWritable value1 = new IntArrayWritable();	
//			value1.set(SimplePolymerChainCodec.encodePolymerChain(polymerType.ordinal(), coords, sequence, gaps));
			long start = System.nanoTime();
			value1.set(SimplePolymerChainCDF53Codec.encodePolymerChain(polymerType.ordinal(), coords, sb.toString(), transform));
		    metrics[1] += System.nanoTime() - start;
			writer.append(key1, value1);
			
			metrics[0] += value1.get().length;
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
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(true);
		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(false);
		cache.setFileParsingParams(params);
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
		//		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		return cache;
	}
}
