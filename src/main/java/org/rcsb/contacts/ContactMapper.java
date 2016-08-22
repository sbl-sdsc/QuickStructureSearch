package org.rcsb.contacts;

import java.util.ArrayList;
import java.util.List;

import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.Grid;
import org.biojava.nbio.structure.contact.GroupContact;
import org.biojava.nbio.structure.contact.GroupContactSet;

import scala.Tuple2;

/**
 * This class generates a tuple of pdbId, pairwise residue contacts between protein chains. The cutoff
 * distance is the minimum distance between any two atoms in a group (residue).
 * @author Peter Rose
 *
 */
public class ContactMapper implements PairFlatMapFunction<Tuple2<String,Structure>,String,String> {
	private static final long serialVersionUID = 4779823229834914377L;
	private double cutoff;

	public ContactMapper(double cutoff) {
		this.cutoff = cutoff;
	}
	
	@Override
	public Iterable<Tuple2<String, String>> call(Tuple2<String, Structure> t) throws Exception {
		List<Tuple2<String,String>> results = new ArrayList<>();

		List<Chain> chains = t._2.getChains();

		for (int i = 0; i < chains.size()-1; i++) {
			Atom[] iAtoms = getChainAtoms(chains.get(i));
			String iChainId = chains.get(i).getName();
			if (iAtoms.length == 0) {
				continue;
			}
			for (int j = i + 1; j < chains.size(); j++) {
				Atom[] jAtoms = getChainAtoms(chains.get(j));
				String jChainId = chains.get(j).getName();
				if (jAtoms.length == 0) {
					continue;
				}
				Grid g = new Grid(cutoff);
				g.addAtoms(iAtoms, jAtoms);
				AtomContactSet contacts = g.getContacts();
				GroupContactSet groupContacts = new GroupContactSet(contacts);
	
				// save results as a comma separated string
				for (GroupContact c: groupContacts) {
					String result = iChainId + "," + 
							c.getPair().getFirst().getResidueNumber() + ","
							+ jChainId + "," + c.getPair().getSecond().getResidueNumber();
					results.add(new Tuple2<String,String>(t._1, result));
					System.out.println(result);
				}
			}
		}
		return results;
	}
	
	private static Atom[] getChainAtoms(Chain chain) {
		Atom[] atoms = new Atom[0];
        if (! chain.getEntityType().equals(EntityType.POLYMER)) {
        	return atoms;
        }
        if (! StructureTools.isProtein(chain)) {
        	return atoms;
        }
        
        return StructureTools.getAllAtomArray(chain);
	}
}
