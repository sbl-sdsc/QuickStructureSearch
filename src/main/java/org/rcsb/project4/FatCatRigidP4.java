package org.rcsb.project4;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.fatcat.calc.AFPCalculator;
import org.biojava.nbio.structure.align.fatcat.calc.AFPChainer;
import org.biojava.nbio.structure.align.fatcat.calc.AFPOptimizer;
import org.biojava.nbio.structure.align.fatcat.calc.AFPPostProcessor;
import org.biojava.nbio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.nbio.structure.align.fatcat.calc.SigEva;
import org.biojava.nbio.structure.align.model.AFP;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;


public class FatCatRigidP4{
	
	public static final String VERSION = "1.1";

	public static String algorithmName = "jFatCat_rigid";

	FatCatParameters params;

	public FatCatRigidP4(){
		super();
		params = new FatCatParameters();
	    params.setMaxTra(0);		
	}

	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param)
	throws StructureException {

		if ( ! (param instanceof FatCatParameters)){
			throw new IllegalArgumentException("FatCat algorithm needs FatCatParameters object as argument.");
		}

		params = (FatCatParameters) param;		
		
		long tstart = System.currentTimeMillis();
		
		AFPChain afpChain = new AFPChain();
	    afpChain.setCa1Length(ca1.length);
	    afpChain.setCa2Length(ca2.length);
		
	    AFPCalculator.extractAFPChains(params, afpChain,ca1, ca2);
	    
	    AFPCalculator.sortAfps(afpChain,ca1,ca2);
	    
	    rChainAfp(params, afpChain,ca1,ca2);
	    
	    long end = System.currentTimeMillis();
		afpChain.setCalculationTime(end-tstart);
		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(VERSION+"");
		return afpChain;
	}
	
	private static void rChainAfp(FatCatParameters params, AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException{
	     params.setMaxTra(0);
	     afpChain.setMaxTra(0);
	     chainAfp(params,afpChain,ca1,ca2);
	}
	
	private static void chainAfp(FatCatParameters params,AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException{
	     
		// we don;t want to rotate input atoms, do we?
		  Atom[] ca2clone = StructureTools.cloneCAArray(ca2);
		  
	     List<AFP> afpSet = afpChain.getAfpSet();

	     int afpNum = afpSet.size();

	     if ( afpNum < 1)
	        return;

	     //run AFP chaining

	     AFPChainer.doChainAfp(params,afpChain ,ca1,ca2);

	     int afpChainLen = afpChain.getAfpChainLen();
	     
	     if(afpChainLen < 1)     {
	        
	        afpChain.setShortAlign(true);
	        return;
	     } //very short alignment
	     
	     // do post processing
	     
	     AFPPostProcessor.postProcess(params, afpChain,ca1,ca2);
	          
	     // Optimize the final alignment 
	     
	     AFPOptimizer.optimizeAln(params, afpChain,ca1,ca2);

	     AFPOptimizer.blockInfo( afpChain);

	     AFPOptimizer.updateScore(params,afpChain);

	     AFPAlignmentDisplay.getAlign(afpChain,ca1,ca2);
	                                            
	     Group[] twistedPDB = AFPTwister.twistPDB(afpChain, ca1, ca2clone);
	     
	     SigEva sig =  new SigEva();
	     double probability = sig.calSigAll(params, afpChain);
	     afpChain.setProbability(probability);
	     double normAlignScore = sig.calNS(params,afpChain);
	     afpChain.setNormAlignScore(normAlignScore);

	     return;

	  }

}
