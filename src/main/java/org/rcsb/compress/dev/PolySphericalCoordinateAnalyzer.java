package org.rcsb.compress.dev;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.vecmath.Point3d;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IOUtils;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.BZip2Codec;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.HTTPConnectionTools;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
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
import org.rcsb.hadoop.io.SimplePolymerChain;
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
public class PolySphericalCoordinateAnalyzer {
	private static final Logger logger = LoggerFactory.getLogger(PolySphericalCoordinateAnalyzer.class);
	private static AtomCache cache = initializeCache();
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";
	

	public static void main(String[] args) throws Exception {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

//		String uri = args[0] 
//				+ "_"
//				+ timeStamp 
//				+ ".seq";
		
		
		List<String> subset = new ArrayList<>(getAll());
//		List<String> pdbIds = subset;
		List<String> pdbIds = subset.subList(90000, 90001);

		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/peter/Data/PDB/");

		long start = System.nanoTime();
		
		String uri = "";
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

		Atom[] atoms = StructureTools.getAllAtomArray(s);
		Map<Atom, Integer> atomMap = new LinkedHashMap<>();

		for (int i = 0; i < atoms.length; i++) {
			for (Bond bond: atoms[i].getBonds()) {
				Atom a = bond.getAtomA();
				Atom b = bond.getAtomB();
				int aIndex = atomMap.get(a);
				int bIndex = atomMap.get(b);
				if (aIndex < bIndex) {
					System.out.println(aIndex + " -"  + bIndex);
				}
				
			}
		}
	
		return 0;
	}
	
	private static int toUnsignedInt(final int n) {
	    return (n << 1) ^ (n >> 31);
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
		params.setCreateAtomBonds(true);
		cache.setFileParsingParams(params);
		DownloadChemCompProvider p = new DownloadChemCompProvider();
		ChemCompGroupFactory.setChemCompProvider(p);
		p.checkDoFirstInstall();
		return cache;
	}
}
