package org.rcsb.ProteinLigandInteractionSearch;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.rcsb.GetRepresentatives;

import com.sangupta.bloomfilter.BloomFilter;
import com.sangupta.bloomfilter.impl.InMemoryBloomFilter;
/**
 * This class creates a hadoop sequence file of protein-ligand interactions,
 * the file contains chain id, sequence number, insertion code, residue name, atom name and element name
 * @author Hinna Shabir
 *
 */

public class CreateBloomFilter {
	private static AtomCache cache = initializeCache();
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";
	/**
	 * 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		StructureIO.setAtomCache(cache);
		long start = System.nanoTime();
		int failure = 0;
		int success = 0;
		int chains = 0;
		
//	    int n = 1661777; // for entire PDB
	    int n = 3000; // per structure
	    double falsePositiveRate = 0.01;
	    
	    Map<String, BloomFilter<String>> bloomMap = new HashMap<>();
	   
//	    BloomFilter<String> filter = new InMemoryBloomFilter<String>(n, falsePositiveRate);


		List<String> subset = new ArrayList<>(GetRepresentatives.getAll());
		//		List<String> pdbIds = subset;
//		List<String> pdbIds = subset.subList(0, 100);
//			List<String> pdbIds = subset.subList(0, 1000);
		List<String> pdbIds = subset.subList(105000, 106000); // 90-95K: 1719, 100-105K: 2095, 105-110K:2879
//		List<String> pdbIds = subset;

//		Set<String> mySet =new HashSet<String>();
		int maxInteraction = 0;
		for (String pdbId: pdbIds) {
			Set<String> mySet =new HashSet<String>();
			System.out.println(pdbId);
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
			List<Group> hetgroups= new ArrayList<Group>();
			List<Group> proteingroups= new ArrayList<Group>();
			List<Chain> chains1 = s.getChains();
			for (Chain c: chains1) {
				List<Group> groups = c.getAtomGroups();
				for (Group group : groups) {
					if (group.getType().equals(GroupType.HETATM) && !group.isWater()) {
						hetgroups.add(group);
					}
					if (group.getType().equals(GroupType.AMINOACID)) {
						proteingroups.add(group);
					}
				}
				chains++;
			}
			Character inscode = null;
			Character inscode2= null;
			for (Group g1: proteingroups){
				for (Atom atom1: g1.getAtoms()) {
					if (!(atom1.getElement().equals(Element.C) ||  atom1.getElement().equals(Element.H)) ){ // Ignore Carbon and Hydrogen atoms in proteins
						Point3d p1= new Point3d (atom1.getCoords());
						for (Group g2 : hetgroups) { 
							for (Atom atom2: g2.getAtoms()) {
								if (!(atom2.getElement().equals(Element.C) ||  atom2.getElement().equals(Element.H))){ // Ignore Carbon and Hydrogen atoms in Ligands
									Point3d p2= new Point3d (atom2.getCoords());
									double d=p1.distance(p2);
									int idist= (int) Math.round(2*d);
									if (idist<=10){ // max distance = 10/2 = 5 A
										ResidueNumber r = g1.getResidueNumber();
										ResidueNumber r2 = g2.getResidueNumber();
										if(r.getInsCode() == null) {
											inscode=' ';
										}
										else{
											inscode=r.getInsCode();
										}
										if(r2.getInsCode() == null) {
											inscode2=' ';
										}
										else{
											inscode2=r.getInsCode();
										}
								//		String label=idist + g1.getChemComp().getId() + g2.getChemComp().getId() + atom1.getName() + atom2.getName();
										String label=idist + g2.getChemComp().getId() + atom1.getName() + atom2.getName();
										mySet.add(label);
	//									System.out.println(label);
										BloomFilter<String> filter = bloomMap.get(g1.getChemComp().getId());
										if (filter == null) {
											filter = new InMemoryBloomFilter<String>(n, falsePositiveRate);
											bloomMap.put(g1.getChemComp().getId(), filter);
										}
									    filter.add(label);
									}
								}
							}
						}
					}     
				}
			}
			maxInteraction = Math.max(mySet.size(), maxInteraction);
		} 

		System.out.println("Total structures: " + pdbIds.size());
	    System.out.println("BloomFilter nbits: " + bloomMap.get("ARG").getNumberOfBits());
	    System.out.println("Maximal number of interactions per structure: " + maxInteraction);
	    System.out.println("BloomFilter contains: 9ARGFMNNEO2P: " + bloomMap.get("ARG").contains("9FMNNEO2P"));
	    System.out.println("BloomFilter contains: 1ARGFMNNEO2P: " + bloomMap.get("ARG").contains("1FMNNEO2P"));
	    System.out.println("BloomFilter contains: 10CYSSF4SGFE2: " + bloomMap.get("CYS").contains("10SF4SGFE2"));
	    System.out.println("BloomFilter contains: 9CYSSF4SGFE3: " + bloomMap.get("CYS").contains("9SF4SGFE3"));
		System.out.println("Success: " + success);
		System.out.println("Failure: " + failure);
		System.out.println("Chains: " + chains);
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

//	/**
//	 * Returns the current list of all PDB IDs.
//	 * @return representatives set of all PDB IDs.
//	 */
//	public static SortedSet<String> getAll() {
//		SortedSet<String> representatives = new TreeSet<String>();
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
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//
//		return representatives;
//	}
	/**
	 * 
	 * @return
	 */
	private static AtomCache initializeCache() {
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
//		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		//params.setParseCAOnly(true);
//		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(false);
		cache.setFileParsingParams(params);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		return cache;
	}
}