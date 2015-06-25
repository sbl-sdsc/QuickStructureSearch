package org.rcsb.project4;

import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Point3d;
import org.apache.spark.Accumulator;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.fatcat.calc.AFPChainer;
import org.biojava.nbio.structure.align.fatcat.calc.AFPOptimizer;
import org.biojava.nbio.structure.align.fatcat.calc.AFPPostProcessor;
import org.biojava.nbio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.nbio.structure.align.fatcat.calc.SigEva;
import org.biojava.nbio.structure.align.model.AFP;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.jama.Matrix;
import org.rcsb.structuralAlignment.SuperPositionQCP;

/**
 * This class calculates the TM score of two protein chains using the FatCatRigid.
 * Use the new QCP method to improve the speed of getting RMSD
 * 
 * @author Chris Li
 */
public class FatCatRigidP4{
	// Parameters of FatCatRigid
	FatCatParameters params;
	// Timers for parallel time counting
	List<Accumulator<Long>> timers;

	public FatCatRigidP4(){
		params = new FatCatParameters();
		params.setMaxTra(0);		
	}

	/**
	 * Align two protein chains
	 * @param ca1
	 * @param ca2
	 * @param timers
	 * @return
	 * @throws StructureException
	 */
	public AFPChain align(Atom[] ca1, Atom[] ca2, List<Accumulator<Long>> timers) throws StructureException {			
		// Load timers
		this.timers = timers;
		AFPChain afpChain = new AFPChain();
		afpChain.setCa1Length(ca1.length);
		afpChain.setCa2Length(ca2.length);

		extractAFPChains(afpChain,ca1, ca2);
		sortAfps(afpChain,ca1,ca2);
		chainAfp(afpChain,ca1,ca2);

		return afpChain;
	}

	/**
	 * Chain AFP, mainly use the doChainAfp (major time consuming)
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @throws StructureException
	 */
	private void chainAfp(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException{
		params.setMaxTra(0);
		afpChain.setMaxTra(0);
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = StructureTools.cloneCAArray(ca2);
		List<AFP> afpSet = afpChain.getAfpSet();
		int afpNum = afpSet.size();
		if (afpNum < 1)
			return;
		//run AFP chaining
		// Timer for doChainAfp
		long startTime = System.nanoTime();
		AFPChainer.doChainAfp(params,afpChain ,ca1,ca2);
		timers.get(1).add(System.nanoTime() - startTime);
		
		int afpChainLen = afpChain.getAfpChainLen();
		if (afpChainLen < 1)     {
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
		AFPTwister.twistPDB(afpChain, ca1, ca2clone);

		SigEva sig =  new SigEva();
		double probability = sig.calSigAll(params, afpChain);
		afpChain.setProbability(probability);
		double normAlignScore = sig.calNS(params,afpChain);
		afpChain.setNormAlignScore(normAlignScore);
		return;
	}

	/**
	 * Extract AFP
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @throws StructureException
	 */
	private void extractAFPChains(AFPChain afpChain,Atom[] ca1,Atom[] ca2) throws StructureException {
		List<AFP> afpSet = new ArrayList<AFP>();
		afpChain.setAfpSet(afpSet);
		int p1, p2;
		double filter1;
		double rmsd = 0;
		Matrix r = new Matrix(3,3);
		Atom   t = new AtomImpl();
		int sparse = params.getSparse();
		int maxTra = params.getMaxTra();
		int fragLen = params.getFragLen();
		double disFilter = params.getDisFilter();
		double rmsdCut = params.getRmsdCut();
		double badRmsd = params.getBadRmsd();
		double fragScore = params.getFragScore();
		int add = sparse + 1; //if add > 1, use sparse sampling
		int minLen = 0;
		int prot1Length = ca1.length;
		int prot2Length = ca2.length;
		if(prot1Length < prot2Length)
			minLen = prot1Length;
		else
			minLen = prot2Length;
		afpChain.setMinLen(minLen);
		afpChain.setBlockResList(new int[maxTra+1][2][minLen]);
		afpChain.setFocusRes1(new int[minLen]);
		afpChain.setFocusRes2(new int[minLen]);

		for(p1 = 0; p1 < prot1Length - fragLen; p1 += add )    {
			for(p2 = 0; p2 < prot2Length - fragLen; p2 += add)     {
				filter1 = getEnd2EndDistance(ca1, ca2, p1, p1 + fragLen - 1, p2, p2 + fragLen - 1);
				//difference between end-to-end distances
				if(filter1 > disFilter) { 
					continue;
				}
				boolean filter2 = filterTerminal(ca1,ca2, p1, p1 + fragLen - 1, p2, p2 + fragLen - 1, fragLen, minLen);
				if(filter2)     {
					continue;
				} //be cautious to use this filter !!

				// Timers for getRmsd
				long startTime = System.nanoTime();
				// here FATCAT does a a jacobi transformation
				//rmsd = kearsay(fragLen, ca1[p1], ca2[p2], r, t);
				// Use the BioJava SVD instead:
				//rmsd = getRmsd(ca1,ca2,fragLen, p1,p2,r,t);
				// Use QCP instead:
				rmsd = getRmsdP3d(ca1,ca2,fragLen, p1,p2);
				timers.get(0).add(System.nanoTime() - startTime);

				if(rmsd < rmsdCut)      {
					AFP afptmp = new AFP();
					afptmp.setP1(p1);
					afptmp.setP2(p2);
					afptmp.setFragLen(fragLen);
					afptmp.setRmsd(rmsd);
					afptmp.setM(r);
					afptmp.setT(t.getCoords());
					afptmp.setScore(scoreAfp(afptmp,badRmsd,fragScore));
					afpSet.add(afptmp);
				}
			}
		}
	}
	
	/**
	 * Sort AFP
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 */
	private void sortAfps(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		List<AFP> afpSet = afpChain.getAfpSet();
		// Get length
		int pro1Len = ca1.length;
		int pro2Len = ca2.length;

		afpChain.setAfpIndex(new int[pro1Len][pro2Len]); //the index of (i,j) pair in AFP list, otherwise -1
		afpChain.setAfpAftIndex(new int[pro1Len][pro2Len]);  //the index of AFP (i,j*) nearest to (i,j), j*<j. if a AFP exits for (i,j), it equals to afpIndex
		afpChain.setAfpBefIndex(new int[pro1Len][pro2Len]); //the index of AFP (i,j*) nearest to (i,j), j*>j. if a AFP exits for (i,j), it equals to afpIndex

		int[][] afpIndex       = afpChain.getAfpIndex();
		int[][] afpAftIndex    = afpChain.getAfpAftIndex();
		int[][] afpBefIndex    = afpChain.getAfpBefIndex();

		for(int i = 0; i < pro1Len; i ++)   {
			for(int j = 0; j < pro2Len; j ++)   {
				afpIndex[i][j] = afpAftIndex[i][j] = afpBefIndex[i][j] = -1;
			}
		}

		//index the AFP for easy extraction of compatible AFPs
		int afpNum = afpSet.size();
		int b0 = 0;
		for(int a = 0; a < afpNum; a ++)    {
			if(a == afpNum - 1 || afpSet.get(a).getP1() != afpSet.get(a+1).getP1())   {
				int i = afpSet.get(a).getP1();
				for(int b = b0; b <= a; b ++)       {
					int j = afpSet.get(b).getP2();
					afpIndex[i][j]=b ;
					afpBefIndex[i][j]=b;
					afpAftIndex[i][j]=b;
					if(afpSet.get(b).getP1() != i)    {
						System.err.println(String.format("Warning: wrong afp index %d %d\n", i, afpSet.get(b).getP1()));
						return;
					}
				}
				for(int k = 1; k < pro2Len; k ++)   {
					if( afpBefIndex[i][k] == -1){
						afpBefIndex[i][k] = afpBefIndex[i][k-1];
					}
				}
				for(int k = pro2Len - 2; k >= 0; k --)      {
					if(afpAftIndex[i][k] == -1) {
						afpAftIndex[i][k] =  afpAftIndex[i][k+1];
					}
				}
				b0 = a + 1;
			}
		}
	}
	
	/**
	 * Get the end to end distance of two chains
	 * @param ca1
	 * @param ca2
	 * @param p1b
	 * @param p1e
	 * @param p2b
	 * @param p2e
	 * @return
	 */
	private double getEnd2EndDistance(Atom[] ca1, Atom[] ca2, int p1b, int p1e, int p2b, int p2e)
	{
		double min = 99;
		double dist1 = Calc.getDistance(ca1[p1b], ca1[p1e]);
		double dist2 = Calc.getDistance(ca2[p2b], ca2[p2e]);
		min = dist1 - dist2;
		return Math.abs(min);
	}
	
	/**
	 * Check if the filter is terminal
	 * @param ca1
	 * @param ca2
	 * @param p1b
	 * @param p1e
	 * @param p2b
	 * @param p2e
	 * @param fragLen
	 * @param minLen
	 * @return
	 */
	private boolean filterTerminal(Atom[] ca1, Atom[] ca2, int p1b, int p1e, int p2b, int p2e, int fragLen, int minLen)
	{
		int d1 = (p1b < p2b)?p1b:p2b;
		int d2 = (ca1.length - p1e) < (ca2.length - p2e)?(ca1.length - p1e):(ca2.length - p2e);
		int d3 = d1 + d2 + fragLen; //maximum alignment length from current AFP
		/// DO NOT DO Math.round() this will give different results to FATCAT....
		int d4 = (int)(0.3 * minLen);
		return d3 < d4;
	}
	
	/**
	 * Use BioJava SVD to get RMSD
	 * @param ca1
	 * @param ca2
	 * @param fragLen
	 * @param p1
	 * @param p2
	 * @param m
	 * @param t
	 * @return
	 * @throws StructureException
	 */
	@SuppressWarnings("unused")
	private double getRmsd(Atom[] ca1, Atom[] ca2, int fragLen, int p1, int p2, Matrix m, Atom t) throws StructureException {
		double rmsd = 99.9;
		Atom[] catmp1 = getFragment(ca1, p1, fragLen,false);
		Atom[] catmp2 = getFragment(ca2, p2, fragLen,true); // clone the atoms for fragment 2 so we can manipulate them...
		if ( catmp1 == null) {
			System.err.println("could not get fragment for ca1 " + p1 + " " + fragLen );
			return rmsd;
		}
		if ( catmp2 == null) {
			System.err.println("could not get fragment for ca2 " + p2 + " " + fragLen );
			return rmsd;
		}
		SVDSuperimposer svd = new SVDSuperimposer(catmp1, catmp2);

		m = svd.getRotation();
		t = svd.getTranslation();

		for (Atom a : catmp2){
			Calc.rotate(a,m);
			Calc.shift(a,t);
		}
		rmsd = SVDSuperimposer.getRMS(catmp1,catmp2);
		return rmsd;
	}
	
	/**
	 * Use the new QCP SuperPosition method to get RMSD
	 * @param ca1
	 * @param ca2
	 * @param fragLen
	 * @param p1
	 * @param p2
	 * @return
	 * @throws StructureException
	 */
	private double getRmsdP3d(Atom[] ca1, Atom[] ca2, int fragLen, int p1, int p2) throws StructureException {
		double rmsd = 99.9;
		Point3d[] catmp1 = getFragmentP3d(ca1, p1, fragLen,false);
		Point3d[] catmp2 = getFragmentP3d(ca2, p2, fragLen,true); // clone the atoms for fragment 2 so we can manipulate them...
		if ( catmp1 == null) {
			System.err.println("could not get fragment for ca1 " + p1 + " " + fragLen );
			return rmsd;
		}
		if ( catmp2 == null) {
			System.err.println("could not get fragment for ca2 " + p2 + " " + fragLen );
			return rmsd;
		}
		SuperPositionQCP qcp = new SuperPositionQCP();
		qcp.set(catmp1, catmp2);
		rmsd = qcp.getRmsd();
		return rmsd;
	}
	
	/** 
	 * Get a fragment of a chain
	 * @param caall
	 * @param pos
	 * @param fragmentLength
	 * @param clone
	 * @return
	 */
	private Atom[] getFragment(Atom[] caall, int pos, int fragmentLength , boolean clone){
		if ( pos+fragmentLength > caall.length)
			return null;
		Atom[] tmp = new Atom[fragmentLength];
		for (int i=0;i< fragmentLength;i++){
			if (clone){
				tmp[i] = (Atom)caall[i+pos].clone();
			} else {
				tmp[i] = caall[i+pos];
			}
		}
		return tmp;
	}
	
	/**
	 * Get a fragment of a chain, and return it as Point3d array
	 * @param caall
	 * @param pos
	 * @param fragmentLength
	 * @param clone
	 * @return
	 */
	private Point3d[] getFragmentP3d(Atom[] caall, int pos, int fragmentLength , boolean clone){
		if ( pos+fragmentLength > caall.length)
			return null;
		Atom[] tmp = new Atom[fragmentLength];
		for (int i=0;i< fragmentLength;i++){
			if (clone){
				tmp[i] = (Atom)caall[i+pos].clone();
			} else {
				tmp[i] = caall[i+pos];
			}
		}
		// Change to Point3d
		Point3d[] chain = new Point3d[tmp.length];
		for (int i=0;i<tmp.length;i++) {
			chain[i] = new Point3d(tmp[i].getX(),tmp[i].getY(),tmp[i].getZ());
		}
		return chain;
	}
	
	/**
	 * Get the score
	 * @param afp
	 * @param badRmsd
	 * @param fragScore
	 * @return
	 */
	private double scoreAfp(AFP afp, double badRmsd, double fragScore)
	{
		//longer AFP with low rmsd is better
		double s, w;
		w = afp.getRmsd() / badRmsd;
		w = w * w;
		s = fragScore * (1.0 - w);
		return s;
	}
}
