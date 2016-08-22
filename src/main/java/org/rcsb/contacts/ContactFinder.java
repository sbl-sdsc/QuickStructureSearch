package org.rcsb.contacts;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.spark.api.java.JavaPairRDD;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.contact.AtomContact;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.Grid;
import org.biojava.nbio.structure.contact.GroupContact;
import org.biojava.nbio.structure.contact.GroupContactSet;
import org.biojava.spark.utils.BiojavaSparkUtils;

public class ContactFinder {

	/**
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		// Starter counter
		Long startTime = System.currentTimeMillis();
		
		// Read file with PDB IDs
		Set<String> pdbIds = getPdbIds(args[0]);
		System.out.println("Dataset size: " + pdbIds.size());
		System.out.println(pdbIds);
				
		double cutoff = 5.0;
		
		BiojavaSparkUtils
				.getBiojavaRdd("/users/peter/Data/mmtf/full")
		        .filter(t -> pdbIds.contains(t._1))
		        .flatMapToPair(new ContactMapper(cutoff))
		        .repartition(1)
		        .saveAsTextFile(args[1]);
	}
	
	private static Set<String> getPdbIds(String fileName) throws FileNotFoundException, IOException {
		Set<String> pdbIds = new TreeSet<>();
		try(BufferedReader br = new BufferedReader(new FileReader(fileName))) {
			String line = br.readLine();
			while (line != null) {
				pdbIds.add(line.trim());
				line = br.readLine();
			}
		}
		return pdbIds;
	}
}
