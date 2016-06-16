package org.rcsb.ProteinLigandInteractionSearch;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashSet;
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
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.HTTPConnectionTools;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
/**
 * This class creates a hadoop sequence file of protein-ligand interactions,
 * the file contains chain id, sequence number, insertion code, residue name, atom name and element name
 * @author Hinna Shabir
 *
 */

public class WriteSeqFile {
	private static AtomCache cache = initializeCache();
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";
/**
 * 
 * @throws IOException 
 */
	public static void main(String[] args) throws IOException {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		String uri = "/Users/peter/Data/PLInteractions/pl" 
				+ "_"
				+ timeStamp 
				+ ".seq";

		Set<String> pdbIds = getAll();
//		StructureIO.setAtomCache(cache);
//		cache.setPath("/Users/hina/DistanceData/");
		long start = System.nanoTime();
		SequenceFile.Writer writer = null;
		int failure = 0;
		int success = 0;
		int chains = 0;

		try {
					writer =SequenceFile.createWriter(new Configuration(),
					SequenceFile.Writer.file(new Path(uri)),
					SequenceFile.Writer.keyClass(Text.class),
					SequenceFile.Writer.valueClass(Text.class),
					SequenceFile.Writer.compression(SequenceFile.CompressionType.BLOCK, new BZip2Codec()));	

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
										int idist= (int) Math.round(d*10);
										if (idist<=50){
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
											String label=r.getChainId() + "," + r2.getChainId() + "," + r.getSeqNum() + "," +  // store the eight selected features as a string for each protein-ligand pair
											r2.getSeqNum() + "," + inscode+ "," + inscode2+ ","+g1.getChemComp().getId()+ "," +
											g2.getChemComp().getId()+ "," + atom1.getName()+ "," + atom2.getName()+ ","
											+ atom1.getElement()+ "," + atom2.getElement()+ ","+ idist;
											mySet.add(label);
										}

									}
								}
							}
						}     
					}
				}

				for (String temp: mySet){	
                    Text key = new Text(temp);
            		Text value = new Text(pdbId);
					writer.append(key, value);

				}
			}	

		} finally {
			IOUtils.closeStream(writer);
		}
		IOUtils.closeStream(writer);
		System.out.println("Total structures: " + pdbIds.size());
		System.out.println("Success: " + success);
		System.out.println("Failure: " + failure);
		System.out.println("Chains: " + chains);
		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}

	/**
	 * Returns the current list of all PDB IDs.
	 * @return representatives set of all PDB IDs.
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
		ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
		return cache;
	}
}