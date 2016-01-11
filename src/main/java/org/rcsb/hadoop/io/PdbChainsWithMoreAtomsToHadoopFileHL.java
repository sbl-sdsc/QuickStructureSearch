package org.rcsb.hadoop.io;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Collection;
import java.util.List;
import java.util.Set;
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

/**
 * This class creates a Hadoop sequence file for protein chains in the PDB. The Hadoop sequence file
 * uses a delta encoding of the PDB coordinates as well as BZIP2 block level compression.
 * 
 * Add CA, N and C atoms, the output file is used for the StructuralAlphabetFingerprint.
 * 
 * @author  Peter Rose, Chris Li
 */
public class PdbChainsWithMoreAtomsToHadoopFileHL {
	private static AtomCache cache = initializeCache();
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";

	public static void main(String[] args) throws Exception {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		String uri = args[0] 
				+ "_"
				+ timeStamp 
				+ ".seq";

		Set<String> pdbIds = getAll();

		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/Chris/Documents/RCSB/Data/Protein_chains/cache/");

		long start = System.nanoTime();
		toSequenceFile(uri, pdbIds, true);
		long end = System.nanoTime();

		System.out.println("Time: " + (end - start)/1E9 + " sec.");
	}

	public static void toSequenceFile(String fileName, Collection<String> pdbIds, boolean verbose)
			throws Exception {

		int failure = 0;
		int success = 0;
		int chains = 0;

		try (SequenceFile.Writer writer = SequenceFile.createWriter(new Configuration(),
				SequenceFile.Writer.file(new Path(fileName)),
				SequenceFile.Writer.keyClass(Text.class),
				SequenceFile.Writer.valueClass(SimplePolymerChain.class),
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

				chains += append(writer, pdbId, s);
			}
			IOUtils.closeStream(writer);
				}

		if (verbose) {
			System.out.println("Total structures: " + pdbIds.size());
			System.out.println("Success: " + success);
			System.out.println("Failure: " + failure);
			System.out.println("Chains: " + chains);
		}
	}

	private static int append(SequenceFile.Writer writer, String pdbId, Structure s)
			throws Exception {

		int chainCount = 0;

		for (Chain c: s.getChains()) {
			List<Group> groups = c.getSeqResGroups(GroupType.AMINOACID);
			System.out.println("seq len: " + c.getSeqResSequence().length());
			System.out.println("seqresgroups.size(): " + c.getSeqResGroups(GroupType.AMINOACID).size());
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

			Point3d[] coords = new Point3d[groups.size()*3];
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
				Atom atomN = null;
				Atom atomC = null;
				if (aminoAcid) {
					atom = ((AminoAcidImpl)g).getCA();
					atomN = ((AminoAcidImpl)g).getN();
					atomC = ((AminoAcidImpl)g).getC();
				} else {
					atom = ((NucleotideImpl)g).getP();
				}
				if (atom == null || atomN == null || atomC == null) {
					gaps++;
				} else {
					coords[i*3] = new Point3d(atom.getCoords());
					coords[i*3+1] = new Point3d(atomN.getCoords());
					coords[i*3+2] = new Point3d(atomC.getCoords());
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

			Text key = new Text(pdbId + "." + c.getChainID());
			SimplePolymerChain value = new SimplePolymerChain();	
			value.setPolymerType(polymerType.ordinal());
			value.setCoordinates(coords);
			value.setSequence(sb.toString());
			writer.append(key, value);
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
		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(false);
		cache.setFileParsingParams(params);
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
		return cache;
	}
}
