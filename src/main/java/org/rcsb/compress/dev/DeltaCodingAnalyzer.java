package org.rcsb.compress.dev;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.rcsb.GetRepresentatives;
import org.rcsb.hadoop.io.SimplePolymerType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class creates a Hadoop sequence file for protein chains in the PDB. The Hadoop sequence file
 * uses a delta encoding of the PDB coordinates as well as BZIP2 block level compression.
 * 
 * 
 * @author  Peter Rose
 */
public class DeltaCodingAnalyzer {
	private static final Logger logger = LoggerFactory.getLogger(DeltaCodingAnalyzer.class);
	private static AtomCache cache = initializeCache();
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";
	

	public static void main(String[] args) throws Exception {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		String uri = args[0] 
				+ "_"
				+ timeStamp 
				+ ".seq";
		
		List<String> subset = new ArrayList<>(GetRepresentatives.getAll());
//		List<String> subset = new ArrayList<>(getAll());
//		List<String> pdbIds = subset;
		List<String> pdbIds = subset.subList(90000, 91000);

		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/peter/Data/PDB/");

		long start = System.nanoTime();
		
		process(uri, pdbIds, true);
	
		long end = System.nanoTime();

		System.out.println("Time: " + (end - start)/1E9 + " sec.");
	}

	public static void process(String fileName,
			Collection<String> pdbIds,
			boolean verbose) throws Exception {

		int failure = 0;
		int success = 0;
		int chains = 0;
		int[] metrics = new int[2];

		DescriptiveStatistics xStats = new DescriptiveStatistics();
		DescriptiveStatistics yStats = new DescriptiveStatistics();
		DescriptiveStatistics zStats = new DescriptiveStatistics();
		DescriptiveStatistics aStats = new DescriptiveStatistics();
			

		for (String pdbId : pdbIds) {
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
			append(s, pdbId, xStats, yStats, zStats, aStats);
		}
		System.out.println("x-values -------------");
		printStatistics(xStats);
		System.out.println("y-values -------------");
		printStatistics(yStats);
		System.out.println("z-values -------------");
		printStatistics(zStats);
		System.out.println("a-values -------------");
		printStatistics(aStats);
		
	}
	
	private static void printStatistics(DescriptiveStatistics stats) {
		System.out.println(stats.toString());
		for (int i = 10; i <= 100; i+=10) {
			System.out.println(stats.getPercentile(i));
		}
	}

	private static int append(Structure s, String pdbId, DescriptiveStatistics xStats, DescriptiveStatistics yStats, DescriptiveStatistics zStats, DescriptiveStatistics aStats) throws Exception {

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
				continue;
			} else if (rna > 0) {
				polymerType = SimplePolymerType.RNA;
				continue;
			} else {
				continue;
			}

	        for (int i = 0, dx = 0, dy = 0, dz = 0; i < coords.size(); i++) {
	        	Point3d p = coords.get(i);
	        	int x = (int) Math.round(p.x*1000);
	        	int y = (int) Math.round(p.y*1000);
	        	int z = (int) Math.round(p.z*1000);
	        	dx = x - dx;
	        	dy = y - dy;
	        	dz = z - dz;
	        	if (i > 0 && Math.abs(dx) < 4000 && Math.abs(dy) < 4000 && Math.abs(dz) < 4000) {
	        			xStats.addValue(dx);
	        			yStats.addValue(dy);
	        			zStats.addValue(dz);	
	        		aStats.addValue(Math.sqrt(dx*dx + dy*dy + dz*dz));
//	        		aStats.addValue(dy);
//	        		aStats.addValue(dz);
	        	}
	 //       	System.out.println("dxyz: " + dx + "," + dy + "," + dz + " ds: " + (dx+dy+dz));
	        	dx = x;
	        	dy = y;
	        	dz = z;

	        }
	
		}
		return chainCount;
	}
	private static int toUnsignedInt(final int n) {
	    return (n << 1) ^ (n >> 31);
	}

	/**
	 * Returns the current list of all PDB IDs.
	 * @return PdbChainKey set of all PDB IDs.
	 */
//	public static SortedSet<String> getAll() {
//		SortedSet<String> representatives = new TreeSet<String>();
//
//		try {
//
//			URL u = new URL(allUrl);
//
//			InputStream stream = HTTPConnectionTools.getInputStream(u, 60000);
//
//			if (stream != null) {
//				BufferedReader reader = new BufferedReader(
//						new InputStreamReader(stream));
//
//				String line = null;
//
//				while ((line = reader.readLine()) != null) {
//					int index = line.lastIndexOf("structureId=");
//					if (index > 0) {
//						representatives.add(line.substring(index + 13, index + 17));
//					}
//				}
//			}
//
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//
//		return representatives;
//	}

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
