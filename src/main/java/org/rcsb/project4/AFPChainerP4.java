package org.rcsb.project4;
/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * http://www.ncbi.nlm.nih.gov/pubmed/14534198
 * </pre>
 * 
 * Thanks to Yuzhen Ye and A. Godzik for granting permission to freely use and redistribute this code.
 *  
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 *
 * Created on Jun 17, 2009
 * Created by Andreas Prlic - RCSB PDB 
 * 
 */
import org.apache.spark.Accumulator;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.nbio.structure.align.model.AFP;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;

import java.util.List;

/** a class to chain AFPs to an alignment
 * 
 * @author Andreas Prlic
 *
 */
public class AFPChainerP4
{
	static List<Accumulator<Long>> classTimers; 
	static double max;
	/**
  // Key function: chain (assembly) the AFPs
  // a AFP (k) is defined as (i, j, k), with i and j are staring points
  // AFP extension (eg. AFP(k-1) -> AFP(k) ) requiring
  // AFP(k-1) < AFP(k)(refer AFP.h definition),
  // ie i(k-1) < i(k) and j(k-1) < j(k)
  // in the figure, only (2) AFP can extend to AFP(k)
  // Key features: a coordination transformation is allowed in the AFP extension
  //                      gap penalties are also considered
  //
  //                                   protein1
  //                 ---------------------------
  //                 |        \                |
  //                 |         \(1)            |
  //                 |     \    \              |
  //                 |      \(2) \             |
  //              p  |       \                 |
  //              r  |   \                     |
  //              o  |    \(3)  \(i,j, k)      |
  //              t  |     \     \             |
  //              e  |            \            |
  //              i  |                         |
  //              n  |                         |
  //              2  ---------------------------
  //                  schematic of AFP chaining
	 * @param timers 
	 */
	public static void doChainAfp(FatCatParameters params, AFPChain afpChain,Atom[] ca1, Atom[] ca2, List<Accumulator<Long>> timers){
		classTimers = timers;
		max = 0;
		List<AFP> afpSet = afpChain.getAfpSet();

		afpChain.setConn(0d);
		afpChain.setDVar(0d);

		int afpNum = afpSet.size();
		if(afpNum <= 0) return;


		int[] twi = afpChain.getTwi();
		if(twi == null) {
			twi = new int[afpNum];
			afpChain.setTwi(twi);
		}

		//transformation, calculated at DoChainAfp, be used in List extraction

		//forward: calculate the score matrix
		boolean isConnected = false;
		int     i, j, j0,  n;
		double  stmp;

		afpChain.setConn(0d);
		afpChain.setDVar(0d);

		double[]  sco = new double[afpNum]; //the score ending at an AFP
		int[] pre = new int[afpNum];    //the previous AFP
		double  maxsco = 0;
		int     maxafp = 0;
		int[] list = new int[afpNum];

		int maxGap = params.getMaxGap();
		int fragLen = params.getFragLen();
		int maxTra = params.getMaxTra();

		Matrix disTable1 = getDisTable(maxGap + 2 * fragLen + 1,ca1);
		Matrix disTable2 = getDisTable(maxGap + 2 * fragLen + 1,ca2);
		

		afpChain.setDisTable1(disTable1);
		afpChain.setDisTable2(disTable2);
		
		for(i = 0; i < afpNum; i ++)    {
			sco[i] = afpSet.get(i).getScore(); //start from itself
			pre[i] = -1;
			twi[i] = 0;
			if ( afpSet.get(i).getP1() < fragLen || afpSet.get(i).getP2() < fragLen)
				n = 0;
			else
				n = getCompatibleAfps(i, list, params, afpChain); //get a compatible list
			//printf("afp %d, compatible %d\n", i, n);
			for(j0 = 0; j0 < n; j0 ++)      {
				j = list[j0];
				if (twi[j] > 0)
					continue;
				isConnected = afpPairConn(j, i, params,afpChain); //note: j, i
				Double conn = afpChain.getConn();
				int t = 0;
				if ( isConnected)
					t=1;
				if(twi[j] + t > maxTra) continue;
				//two many transformation are disfavored
				stmp = sco[j] + afpSet.get(i).getScore() + conn;
				if(stmp > sco[i])       { //considered all previous compatible AFPs
					sco[i] = stmp;
					twi[i] = twi[j] + t;
					pre[i] = j;
				}
			}

			if(maxsco < sco[i])     {
				maxsco = sco[i];
				maxafp = i;
			}
		}
		
		
		int     currafp = maxafp;
		
		//trace-back from maxafp (maxsco)

		afpChain.setAlignScore(maxsco);
		afpChain.setAlignScoreUpdate(maxsco);
		afpChain.setAfpChainTwiNum(0);
	
		traceBack(pre, currafp, twi[currafp],params,afpChain,ca1,ca2);
	}

	private static Matrix getDisTable(int maxlen, Atom[]ca)

	{
		int length = ca.length;
		Matrix dis = new Matrix(length,length);

		int     i, j;
		for(i = 0; i < length; i ++)    {
			dis.set(i,i,0);
			for(j = i + 1;( j < length) && (j <= i + maxlen); j ++)     {
				double val = Calc.getDistance(ca[i],ca[j]);
				dis.set(i,j,val);
				dis.set(j,i,val);
			}
		}
		return dis;

	}

	/*

  derive the compabitle AFP lists for AFP-chaining
  this is important for speeding up the process
  for a given AFP(i1,j1), there are three regions that could be the starting
  point for the compabitle AFPs of AFP(i1,j1)
  //                 a1        a2   a3
  //               i1-G    i1-f-c i1-f+1 i1
  //                 |          |   |   |
  //              ----------------------------
  //              | B               |   |
  //b1  j1-G  -|  ---------------|   |
  //              |  |          |   |   |
  //              |  |     C    | 3 |   |
  //              |  |          |   |   |
  //b2 j1-f-c -|  |--------------|   |
  //              |  |     2    | 1 |   |
  //b3 j1-f+1 -|------------------   |
  //              |                   A |
  //          j1 -|---------------------\
  //              |                      \ (AFP(i1,j1))
  //              -----------------------------
  //
  f: the length of AFPs (we use segments of same length)
  G: g + f, where g is the maximum allowed gaps
  c: the maximum allowed cross-over in AFP-connection,
     here we use c = f, and j1-f-c = j1-2f
  incompatible region A: its AFPs overlap with given AFP(i1,j1)
  incompatible region B: the gaps between its AFP with given AFP is larger than g
  incompatible region C: its AFPs connect with given AFP but cross a given threshold.
  compatible region 1: [i1-f-c,i1-f+1>,[j1-f-c,j1-f+1> or [a2,a3],[b2,b3]
  compatible region 2: [i1-G,i1-f-c],[j1-f-c,j1-f] or [a1,a2],[b2,b3]
  combine 1 and 2    : [i1-G,i1-f],[j1-f-c,j1-f]   or [a1,a3],[b2,b3]
  compatible region 3: [i1-f-c,i1-f],[j1-G,j1-f-c] or [a2,a3],[b1,b2]
  c->misCut
  f->fragLen
  G->fragLen+maxGap->maxGapFrag
	 *
	 *
	 */
	private  static int getCompatibleAfps(int afp, int[] list, FatCatParameters params, AFPChain afpChain){

		int     i, j, i1, j1, f, G, c, a1, a2, a3, b1, b2, b3, s1, s2;

		int fragLen = params.getFragLen();
		int maxGapFrag = params.getMaxGapFrag();
		int misCut = params.getMisCut();
		int maxTra = params.getMaxTra();
		List<AFP> afpSet = afpChain.getAfpSet();

		f = fragLen;
		G = maxGapFrag;
		c = misCut;

		i1 = afpSet.get(afp).getP1();
		j1 = afpSet.get(afp).getP2();
		a3 = i1 - f;
		a2 = a3 - c;
		a1 = i1 - G;
		a2 = a2>0?a2:0;
		a1 = a1>0?a1:0;

		b3 = j1 - f;
		b2 = b3 - c;
		b1 = j1 - G;
		b2 = (b2 > 0)?b2:0;
		b1 = (b1 > 0)?b1:0;

		int[][] afpAftIndex = afpChain.getAfpAftIndex();
		int[][] afpBefIndex = afpChain.getAfpBefIndex();
		int[] twi = afpChain.getTwi();



		int     n = 0;
		//compatible region 1-2, [a1,a3][b2,b3]
				for(i = a1; i <= a3; i ++)      {//i <= a3 instead of i < a3
					s1 = afpAftIndex[i][b2]; //note afpAftIndex, not afpIndex
					if(s1 < 0)      continue;//no AFP for the given i with j > b2
					s2 = afpBefIndex[i][b3]; //afps is sorted by j given a i,it's sparse matrix
					if(s2 < 0)      continue;//no AFP for the given i with j < b3
					for(j = s1; j <= s2; j ++)      { //j <= s2 instead of j < s2
						if(twi[j] <= maxTra)    {
							list[n ++] = j;
						}
					}
				}

				//compatible region 3  [a2,a3][b1,b2]
				for(i = a2; i <= a3; i ++)      {
					s1 = afpAftIndex[i][b1];
					if(s1 < 0)      continue;
					s2 = afpBefIndex[i][b2]; //afps is sorted by j given a i
					if(s2 < 0)      continue;
					//note j < s2, as the cases of j == s2 is alread considered in previous region
					for(j = s1; j < s2; j ++)       {
						if(twi[j] <= maxTra)    {
							list[n ++] = j;
						}
					}
				}

				return n;

	}


	/**
  //Key function: calculate the connectivity of AFP pairs
  //no compatibility criteria is executed
  //note: afp1 is previous to afp2 in terms of the position
     //this module must be optimized
	 *
	 * @param afp1
	 * @param afp2
	 * @return flag if they are connected
	 */
	public static boolean afpPairConn(int afp1, int afp2,  FatCatParameters params, AFPChain afpChain)

	{
		Double conn = afpChain.getConn();
		Double dvar = afpChain.getDVar();

		double misScore = params.getMisScore();
		double maxPenalty = params.getMaxPenalty();
		double disCut = params.getDisCut();
		double gapExtend = params.getGapExtend();
		double torsionPenalty = params.getTorsionPenalty();
		double disSmooth = params.getDisSmooth();

		List<AFP> afpSet = afpChain.getAfpSet();

		int     m = calcGap(afpSet.get(afp2),afpSet.get(afp1));
		int     g = calcMismatch(afpSet.get(afp2),afpSet.get(afp1));


		double  gp = misScore * m;      //on average, penalty for a mismatch is misScore, no modification on score
		if(g > 0)       {
			gp += gapExtend * g;
		}
		if(gp < maxPenalty)     gp = maxPenalty; //penalty cut-off
		//note: use < (smaller) instead of >, because maxPenalty is a negative number

		double  d;
		//TODO
		double d2 = calAfpDisPo2(afp1, afp2,params, afpChain);
		if (d2 < 6)
			d = calAfpDis(afp1, afp2,params, afpChain);
		else 
			d = disCut;
//		d = calAfpDis(afp1, afp2,params, afpChain);
		
		//note: the 'dis' value is numerically equivalent to the 'rms' with exceptions

		boolean     ch = false;
		double  tp = 0.0;
		if(d >= disCut) {
			tp = torsionPenalty;
			ch = true;
		} //use the variation of the distances between AFPs
		else  if(d > disCut - disSmooth)        {
			double  wt = Math.sqrt((d - disCut + disSmooth) / disSmooth);
			//using sqrt: penalty increase with dis more quicker than linear function
			tp = torsionPenalty * wt;
		}

		dvar = d;
		conn = tp + gp;

		afpChain.setConn(conn);
		afpChain.setDVar(dvar);
		return ch;
	}

	/**
	 * return the gaps between this and afp
	 * requiring afp1 >  afp2
	 * ( operator % in C)
	 */

	private static int  calcGap(AFP afp1 , AFP afp2)
	{
		int     g = (afp1.getP1() - afp2.getP1()) - (afp1.getP2() - afp2.getP2());
		if(g < 0)       g = -g;
		return g;
	}

	/**
	 * return the mis-matched between this and afp
	 *     requiring  AFP1 > afp2
	 *     (operator / in C)
	 * @param afp1
	 * @param afp2
	 * @return
	 */

	//--------------------------------------------
	private static int calcMismatch(AFP afp1, AFP afp2)
	{
		int     l1 = afp1.getP1() - afp2.getP1() - afp2.getFragLen();
		int     l2 = afp1.getP2() - afp2.getP2() - afp2.getFragLen();
		return (l1 > l2?l2:l1);
	}

	/**
  //return the root mean square of the distance matrix between the residues
  //from the segments that form the given AFP list
  //this value can be a measurement (2) for the connectivity of the AFPs
  //and its calculation is quicker than the measurement (1), rmsd
  //currently only deal with AFP pair
//
//           |-d1--|
//          |--d2---|
//         |---d3----|
  //-----------------------------------------------------------------------
  //this module is optimized
	 *
	 * @param afp1
	 * @param afp2
	 * @return
	 */
	private static double calAfpDis(int afp1, int afp2, FatCatParameters params, AFPChain afpChain)
	{

		List<AFP> afpSet = afpChain.getAfpSet();

		Matrix disTable1 = afpChain.getDisTable1();
		Matrix disTable2 = afpChain.getDisTable2();

		int fragLen = params.getFragLen();
		double afpDisCut = params.getAfpDisCut();
		double disCut = params.getDisCut();
		double fragLenSq = params.getFragLenSq();

		int     i, j, ai, bi, aj, bj;
		double  d;
		double  rms = 0;
		for(i = 0; i < fragLen; i ++)   {
			ai = afpSet.get(afp1).getP1() + i;
			bi = afpSet.get(afp1).getP2() + i;
			for(j = 0; j < fragLen; j ++)   {
				aj = afpSet.get(afp2).getP1() + j;
				bj = afpSet.get(afp2).getP2() + j;
				d = disTable1.get(aj,ai) - disTable2.get(bj,bi);
				rms += d * d;
				if(rms > afpDisCut)     { return (disCut); }
			}
		}
		return (Math.sqrt(rms / fragLenSq));
	}
	
	@SuppressWarnings("unused")
	private static double calAfpDisPo(int afp1, int afp2, FatCatParameters params, AFPChain afpChain)
	{

		List<AFP> afpSet = afpChain.getAfpSet();

		Matrix disTable1 = afpChain.getDisTable1();
		Matrix disTable2 = afpChain.getDisTable2();

		int fragLen = params.getFragLen();

		int     ai1, bi1, aj1, bj1;
		int 	ai2, bi2, aj2, bj2;
		double  d1, d2, d3, d4;
		double  rms = 0;
		ai1 = afpSet.get(afp1).getP1();
		bi1 = afpSet.get(afp1).getP2();
		aj1 = afpSet.get(afp2).getP1();
		bj1 = afpSet.get(afp2).getP2();
		d1 = disTable1.get(aj1,ai1) - disTable2.get(bj1,bi1);
		ai2 = ai1 + fragLen - 1;
		bi2 = bi1 + fragLen - 1;
		aj2 = aj1 + fragLen - 1;
		bj2 = bj1 + fragLen - 1;
		d2 = disTable1.get(aj2,ai2) - disTable2.get(bj2,bi2);
		d3 = disTable1.get(aj1,ai2) - disTable2.get(bj1,bi2);
		d4 = disTable1.get(aj2,ai1) - disTable2.get(bj2,bi1);
		
		rms += d1 * d1 + d2 * d2 + d3 * d3 + d4 * d4;
		return (Math.sqrt(rms / 4));
	}
	
	private static double calAfpDisPo2(int afp1, int afp2, FatCatParameters params, AFPChain afpChain)
	{

		List<AFP> afpSet = afpChain.getAfpSet();

		Matrix disTable1 = afpChain.getDisTable1();
		Matrix disTable2 = afpChain.getDisTable2();

		int fragLen = params.getFragLen();

		int     ai1, bi1, aj1, bj1;
		int 	ai2, bi2, aj2, bj2, ai3, bi3, aj3, bj3;
		double  d1, d2, d3, d4, d5, d6, d7, d8;
		double  rms = 0;
		ai1 = afpSet.get(afp1).getP1();
		bi1 = afpSet.get(afp1).getP2();
		aj1 = afpSet.get(afp2).getP1();
		bj1 = afpSet.get(afp2).getP2();
		
		ai2 = ai1 + fragLen - 1;
		bi2 = bi1 + fragLen - 1;
		aj2 = aj1 + fragLen - 1;
		bj2 = bj1 + fragLen - 1;
		
		ai3 = ai1 + fragLen/2;
		bi3 = bi1 + fragLen/2;
		aj3 = aj1 + fragLen/2;
		bj3 = bj1 + fragLen/2;
		
		d1 = disTable1.get(aj1,ai1) - disTable2.get(bj1,bi1);
		d2 = disTable1.get(aj2,ai2) - disTable2.get(bj2,bi2);
		d3 = disTable1.get(aj1,ai2) - disTable2.get(bj1,bi2);
		d4 = disTable1.get(aj2,ai1) - disTable2.get(bj2,bi1);
		d5 = disTable1.get(aj3,ai1) - disTable2.get(bj3,bi1);
		d6 = disTable1.get(aj3,ai2) - disTable2.get(bj3,bi2);
		d7 = disTable1.get(aj1,ai3) - disTable2.get(bj1,bi3);
		d8 = disTable1.get(aj2,ai3) - disTable2.get(bj2,bi3);
		
		
		rms += d1 * d1 + d2 * d2 + d3 * d3 + d4 * d4 + d5 * d5 + d6 * d6 + d7 * d7 + d8 * d8;
		return (Math.sqrt(rms / 8));
	}


	/**
	 * derive the optimal chaining of AFPs by trace-back
	 */
	private static void traceBack(int[] pre, int currafp0, int twist, FatCatParameters params, AFPChain afpChain,Atom[] ca1, Atom[]ca2)
	{

		afpChain.setConn(0d);
		afpChain.setDVar(0d);

		int minLen = afpChain.getMinLen();
		List<AFP> afpSet = afpChain.getAfpSet();

		int afpChainLen = 0;

		//trace-back from currafp (maxsco)
		int[]     afpchain    = new int[minLen];
		int[]     afptwibin   = new int[minLen];
		double[]  afptwilist  = new double[minLen];
		int       currafp = currafp0;
		int       s = 0;

		afptwibin[s] = 0;
		afpchain[s ++] = currafp;

		boolean isConnected = false;
		int     prevafp;
		// int     alnlen = afpSet.get(afpchain[s]).getFragLen();

		//Double  conn = new Double(0) ;
		//Double dvar =  new Double(0);

		while((prevafp = pre[currafp]) != -1)   {

			isConnected = afpPairConn(prevafp, currafp, params,afpChain);

			if ( isConnected )
				afptwibin[s - 1] = 1;
			else
				afptwibin[s - 1] = 0;
			Double dvar = afpChain.getDVar();
			afptwilist[s - 1] = dvar;
			//note s - 1: the transformation of prevafp-currafp is recorded in currafp
			currafp = prevafp;
			// alnlen += afpSet.get(currafp).getFragLen();
			afpchain[s ++] = currafp;
		}

		afpChainLen = s;
		afpChain.setAfpChainLen(afpChainLen);

		//first afp without transformation
		if ( isConnected)
			afptwibin[s - 1] = 1;
		else
			afptwibin[s - 1] = 0;

		int[] afpChainList = afpChain.getAfpChainList();
		double[] afpChainTwiBin = afpChain.getAfpChainTwiBin();
		double[] afpChainTwiList = afpChain.getAfpChainTwiList();

		if(afpChainList == null)        {
			afpChainList   = new int[s];
			afpChain.setAfpChainList(afpChainList);
			afpChainTwiBin     = new double[s];
			afpChain.setAfpChainTwiBin(afpChainTwiBin);
			afpChainTwiList = new double[s];
			afpChain.setAfpChainTwiList(afpChainTwiList);
		}

		int afpChainTwiNum = afpChain.getAfpChainTwiNum();

		int     i;
		for(i = 0; i < s; i ++) {
			afpChainList[i]     = afpchain[s - 1 - i];
			afpChainTwiBin[i]   = afptwibin[s - 1 - i];
			afpChainTwiList[i]  = afptwilist[s - 1 - i];
			afpChainTwiNum      += afptwibin[s - 1 - i];
		}

		if(afpChainTwiNum != twist)     {
			System.err.println(String.format("AFPChainer Warning: the twists number is not consistent %d %d\n", afpChainTwiNum, twist));
		}

		double alignScore = afpChain.getAlignScore();

		double  checkscore = afpSet.get(afpChainList[0]).getScore();
		for(i = 1; i < afpChainLen; i ++)       {
			isConnected = afpPairConn(afpChainList[i - 1], afpChainList[i], params,afpChain);
			checkscore = checkscore + afpSet.get(afpChainList[i]).getScore() + afpChain.getConn();
		}
		if(Math.abs(checkscore - alignScore) > 1e-4)        {
			System.err.println(String.format("AFPChainer Warning: fail in alignment score checking %.4f %.4f\n", alignScore, checkscore));
		}

		double  rmsd = calAfpRmsd(afpChainLen, afpChainList,0, afpChain,ca1,ca2);
		afpChain.setChainRmsd(rmsd);
		int     b1 = 0;
		int     bk = 0;
		int     a, b;
		afpChain.setChainLen( 0);
		int chainLen       = afpChain.getChainLen();
		int block2Afp[]    = afpChain.getBlock2Afp();


		double[] blockRmsd = afpChain.getBlockRmsd();
		int[] blockSize = afpChain.getBlockSize();

		block2Afp[0] = 0;
		for(i = 0; i < afpChainLen; i ++)       {
			a = afpChainList[i];
			chainLen += afpSet.get(a).getFragLen();
			if(i > 0)       {
				b = afpChainList[i - 1];
				int misLen = afpChain.getMisLen();
				misLen += calcMismatch(afpSet.get(a),afpSet.get(b));
				afpChain.setMisLen(misLen);
				int gapLen = afpChain.getGapLen();
				gapLen += calcGap(afpSet.get(a),afpSet.get(b));
				afpChain.setGapLen(gapLen);

			}

			if(afpChainTwiBin[i] == 1)   {
				rmsd = calAfpRmsd(i - b1, afpChainList,b1, afpChain, ca1, ca2);

				blockRmsd[bk] = rmsd;
				blockSize[bk] = i - b1;
				b1 = i;

				//System.out.println("block2Afp.length:"+ block2Afp.length + " " + bk + " " + i + " " + afpChain.getMaxTra() );
				block2Afp[++bk] = i; //next-block
			}
		}
		afpChain.setBlock2Afp(block2Afp);
		afpChain.setChainLen(chainLen);

		rmsd = calAfpRmsd(i - b1, afpChainList, b1, afpChain,ca1,ca2);

		// argh this is the block RMSD, not the chain RMSD!
		//afpChain.setChainRmsd(rmsd);
		//rmsd = calAfpRmsd(i - b1, afpChainList[b1]);
		blockSize[bk] = i - b1;
		blockRmsd[bk] = rmsd;

		afpChain.setBlockSize(blockSize);
		afpChain.setBlockRmsd(blockRmsd);
		int blockNum = afpChain.getBlockNum();
		blockNum = ++bk;
		afpChain.setBlockNum(blockNum);
		afpChain.setAfpChainList(afpChainList);
		afpChain.setAfpChainTwiList(afpChainTwiList);
	}

	/**
  //return the rmsd of the residues from the segments that form the given AFP list
  //this value can be a measurement (1) for the connectivity of the AFPs
	 *
	 * @param afpn
	 * @param afpPositions the positions of AFPs to work on. 
	 * @param listStart the starting position in the list of AFPs
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return rmsd
	 */

	protected static double calAfpRmsd(int afpn, int[] afpPositions, int listStart,  AFPChain afpChain,Atom[] ca1, Atom[] ca2)
	{
		int focusResn = AFPTwister.afp2Res(afpChain,afpn, afpPositions, listStart);
		int[] focusRes1 = afpChain.getFocusRes1();
		int[] focusRes2 = afpChain.getFocusRes2();

		double rmsd = getRmsd(focusResn, focusRes1, focusRes2 , afpChain,ca1,  ca2);
		return rmsd;
	}

	/** the the RMSD for the residues defined in the two arrays
	 *
	 * @param focusResn
	 * @param focusRes1
	 * @param focusRes2
	 * @return
	 */
	private static double getRmsd(int focusResn, int[] focusRes1, int[] focusRes2, AFPChain afpChain, Atom[] ca1, Atom[] ca2){
		Atom[] tmp1 = new Atom[focusResn];
		Atom[] tmp2 = new Atom[focusResn];

		for ( int i =0 ; i< focusResn;i++){
			tmp1[i] =       ca1[focusRes1[i]];
			tmp2[i] = (Atom)ca2[focusRes2[i]].clone();
			if (tmp1[i].getCoords() == null){
				System.err.println("tmp1 got null: " +i + " pos: " + focusRes1[i]);
			}
			if (tmp2[i].getCoords() == null){
				System.err.println("tmp1 got null: " +i + " pos: " + focusRes2[i]);
			}
			//XX
			//tmp2[i].setParent((Group) ca2[focusRes2[i]].getParent().clone());
		}
		double rmsd = 99;
		try {
			rmsd = getRmsd(tmp1,tmp2);

		} catch (Exception e){
			e.printStackTrace();
		}
		return rmsd;
	}
	
	/** Calculate the RMSD for two sets of atoms. Rotates the 2nd atom set so make sure this does not cause problems later
	 *
	 *
	 * @param catmp1
	 * @return
	 */
	private static double getRmsd(Atom[] catmp1, Atom[] catmp2) throws StructureException{
		SVDSuperimposer svd = new SVDSuperimposer(catmp1, catmp2);
		Matrix m = svd.getRotation();
		Atom t = svd.getTranslation();

		for (Atom a : catmp2){
			Calc.rotate(a,m);
			Calc.shift(a,t);

		}
		double rmsd = SVDSuperimposer.getRMS(catmp1,catmp2);
		return rmsd;
	}
}
