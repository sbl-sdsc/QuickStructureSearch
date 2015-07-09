package org.rcsb.ProteinLigandInteractionSearch;

import java.io.IOException;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;

public class Ligands {

	public static void main(String[] args) throws IOException, StructureException {
		// TODO Auto-generated method stub
		
        Structure structure = StructureIO.getStructure("1STP");     
        System.out.println(structure);

        List<Chain> chains = structure.getChains();
        System.out.println(chains);
        System.out.println(" # chains: " + chains.size());

        for (Chain c : chains) {

            System.out.println("   Chain: " + c.getChainID() + " # groups with atoms: " + c.getAtomGroups().size());

          //  for (Group g: c.getAtomGroups()){

               // if ( g.getPDBName().equalsIgnoreCase("HEM")) {

            //        System.out.println("   " + g);

                  //  for (Atom a: g.getAtoms()) {

                    //    System.out.println("    " + a);

                    //}
                //}
           // }
            
            
 /*           List<Group> allgroups = c.getAtomGroups();
            for (Group group : allgroups) {
                if ( group instanceof AminoAcid) {
                    AminoAcid aa = (AminoAcid) group;
                    System.out.println(aa.getSecStruc());
                }
            }*/
            
            
        
            //List<Group> groupL = c.getAtomLigands();
            //System.out.println("length: "+ c.getSeqResLength());
            List<Group> groups = c.getAtomGroups();
            for (Group group : groups) {
            	System.out.println("group: " + group);
            //System.out.println("atoms:"+ group.getAtoms());
              if (group.getType().equals(GroupType.HETATM) && !group.isWater()) {
            	//System.out.println("Ligand: " + group); 
            	
            	System.out.println("type: "+ group.getType()); 
            	System.out.println(group.getAtoms());
            	for (Atom a: group.getAtoms()) {
                     if (!(a.getElement().equals(Element.C) ||  a.getElement().equals(Element.H))){
                       //System.out.println("    " + a);
                       //Atom ligandAtms= a;
                       
                       Point3d p1= new Point3d (a.getCoords());
                       
                       for (Group group2 : groups) {
                      	 if (group2.getType().equals(GroupType.AMINOACID) && !group2.isWater()) {
                         	//System.out.println("Protein: " + group); 
                         	System.out.println("Type: "+ group2.getType()); 
                         	for (Atom b: group2.getAtoms()) {
                                  if (!(b.getElement().equals(Element.C) ||  b.getElement().equals(Element.H))){
                                   // System.out.println("    " + b);
                                    //Atom proteinAtms= b;
                                   
                                     Point3d p2= new Point3d (b.getCoords());
                         	          double d=p1.distance(p2);
                         	          int idist= (int) Math.round(d*10);
                         	          //System.out.println(group2.getChemComp().getId()+"-"+ group.getChemComp().getId()+"-"+ b.getElement()+ "-"+ a.getElement()+"-"+ d);
                                      String label=group2.getChemComp().getId()+"-"+ group.getChemComp().getId()+"-"+ b.getName()+ "-"+ a.getName()+"-"+ idist;
                                      System.out.println(label);
                                  }
                                 }
                         	}
                      }        
                     }
            	}   
            }
            }
  	        //System.out.println(d);  
        }
	}
}
