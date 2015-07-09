package org.rcsb.ProteinLigandInteractionSearch;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.vecmath.Point3d;

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
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.HTTPConnectionTools;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.AllChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.rcsb.structuralSimilarity.IntArrayWritable;


public class ProteinLigandDistanceFile {
	private static AtomCache cache = initializeCache();
	private static final int SCALE = 1000; // Factor to convert PDB coordinates to integer values. Do not change this value!
	private static String allUrl = "http://www.rcsb.org/pdb/rest/getCurrent/";

	public static void main(String[] args) throws IOException {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		String uri = "/Users/hina/DistanceData/protein_ligand_distances" 
				+ "_"
				+ timeStamp 
				+ ".seq";

		Set<String> pdbIds = getAll();
		List <String> subset= new ArrayList <>(pdbIds);
		subset = subset.subList(1, 100);
		subset.clear();
		subset.add("1STP");
		StructureIO.setAtomCache(cache);
		cache.setPath("/Users/hina/DistanceData/");

		long start = System.nanoTime();

		Text key = new Text();
		Text value = new Text();	
		SequenceFile.Writer writer = null;

		int failure = 0;
		int success = 0;
		int chains = 0;

		try {
			writer =SequenceFile.createWriter(new Configuration(),
					SequenceFile.Writer.file(new Path(uri)),
					SequenceFile.Writer.keyClass(key.getClass()),
					SequenceFile.Writer.valueClass(value.getClass()),
					SequenceFile.Writer.compression(SequenceFile.CompressionType.BLOCK, new BZip2Codec()));	

			for (String pdbId: subset) {
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
				
				for (Chain c: s.getChains()) {
					System.out.println("   Chain: " + c.getChainID() + " # groups with atoms: " + c.getAtomGroups().size());
					List<Group> groups = c.getAtomGroups();
					for (Group group : groups) {
						System.out.println(group);
						System.out.println("atoms: "+group.getAtoms());
					if (group.getType().equals(GroupType.HETATM) ) {
							//System.out.println("Ligand: " + group); 
							
							System.out.println("type: "+ group.getType()); 
							
							for (Atom aa: group.getAtoms()) {
								System.out.println(" Atoms:  " + aa);
								if (!(aa.getElement().equals(Element.C) ||  aa.getElement().equals(Element.H))){
									System.out.println("    " + aa);
									for (Group group2 : groups) {
										if (group2.getType().equals(GroupType.AMINOACID)) {
											//System.out.println("Protein: " + group); 
											System.out.println("Type: "+ group2.getType()); 
											//System.out.println("atoms: "+group2.getAtoms());
											for (Atom b: group2.getAtoms()) {
												 System.out.println("    " + b);
												if (!(b.getElement().equals(Element.C) ||  b.getElement().equals(Element.H))){
													// System.out.println("    " + b);
													Point3d p1= new Point3d (aa.getCoords());
													Point3d p2= new Point3d (b.getCoords());
													double d=p1.distance(p2);
													int idist= (int) Math.round(d*10);
													String label=group2.getChemComp().getId()+"-"+ group.getChemComp().getId()+"-"+ b.getName()+ "-"+ aa.getName()+"-"+ idist;
													System.out.println(label);   
													key.set(label); 
													value.set(pdbId); 
													writer.append(key, value);
												}
											}
										}
									}     
								}
							}   
						}
				}
					chains++;
				}
			}
		} finally {
			IOUtils.closeStream(writer);
		}

		System.out.println("Total structures: " + pdbIds.size());
		System.out.println("Success: " + success);
		System.out.println("Failure: " + failure);
		System.out.println("Chains: " + chains);
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
