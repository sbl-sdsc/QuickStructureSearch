package org.rcsb.project3;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.apache.spark.broadcast.Broadcast;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.rcsb.structuralAlignment.SuperPositionQCP;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the longest local common subsequence over the length of the chains
 * using SmithWaterman algorithm and Gotoh's improvement
 * 
 * @author Chris Li
 */
public class SmithWatermanWithGeoComp implements AlignmentAlgorithmInterface {

	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data = null;
	private Broadcast<List<Tuple2<String,Point3d[]>>> coords = null;
    /* With different open and extend penalty, this class could function the same as LCS or SmithWaterman
     * LCS: open = extend = 0;
     * SmithWaterman = open = extend = 1;
     */
    // open gap penalty
    private double open = 5;
    // extend gap penalty
    private double extend = 0.5;
    // number of rotate and transform time
    private int rotateTime = 4;

    public SmithWatermanWithGeoComp() {
	}
    
    public SmithWatermanWithGeoComp(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
    
    /***
     * Constructor with setting options
     * @param data
     * @param open
     * @param extend
     */
    public SmithWatermanWithGeoComp(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data, double open, double extend) {
		this.data = data;
		this.open = open;
		this.extend = extend;
	}
    
	@Override
	public void setSequence(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
	
	@Override
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) throws FileNotFoundException {
		Tuple2<String,SequenceFeatureInterface<?>> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,SequenceFeatureInterface<?>> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
			
		SequenceFeatureInterface<?> v1 = t1._2;
		SequenceFeatureInterface<?> v2 = t2._2;
		
		Alignment<?> SWAlignment = getAlignment(v1, v2, open, extend);
		
		// TODO for OneAgainstAll
		Point3d[] c1 = this.coords.getValue().get(tuple._1)._2;
		Point3d[] c2p = this.coords.getValue().get(tuple._2)._2;		
//		Point3d[] c1 = v1.getCoords();
//		Point3d[] c2 = v2.getCoords();
		Point3d[] c2 = new Point3d[c2p.length];
		for (int i = 0; i < c2p.length; i++)
			c2[i] = c2p[i];
//		if (t2._1.equals("1B0B.A")) {
//			PrintWriter writer = new PrintWriter("/Users/Chris/Documents/RCSB/Data/pdb/4HHB.pdb");
//			Atom[] atoms = getCAAtoms(c1);
//			for (Atom a: atoms) {
//				writer.print(toPdb(a).replaceAll("nullGLU", " GLU"));
//			}
//			writer.flush();
//			writer.close();
//			
//			PrintWriter writer2 = new PrintWriter("/Users/Chris/Documents/RCSB/Data/pdb/1B0B.pdb");
//			Atom[] atoms2 = getCAAtoms(c2p);
//			for (Atom a: atoms2) {
//				writer2.print(toPdb(a).replaceAll("nullGLU", " GLU"));
//			}
//			writer2.flush();
//			writer2.close();
//		}
				
		float value = 0;

		// rotate and transform to overlap two chains
		for (int trial = 0; trial < rotateTime; trial++) {
			SWAlignment = rotateAlign(SWAlignment, c1, c2, t1._1,t2._1, trial);
		}

		// FatCat score calculation
		if (SWAlignment != null) {
			Integer[] v1Order = SWAlignment.getSequence1();
			Integer[] v2Order = SWAlignment.getSequence2();
			
			List<Point3d> lp1 = new ArrayList<Point3d>();
			List<Point3d> lp2 = new ArrayList<Point3d>();
			
			for (int i = 0; i < v1Order.length; i++) {
				if (v1Order[i] != null && v2Order[i] != null) {
					lp1.add(c1[v1Order[i]]);
					lp2.add(c2[v2Order[i]]);
				}
			}
	
			int Laln = lp1.size();
			int Lmin = Math.min(c1.length, c2.length);

			double d0 = 1.24 * Math.cbrt(Lmin - 15.) - 1.8;
			double d0sq = d0 * d0;
	
			double sum = 0;
			for (int i = 0; i < Laln; i++) {
				double d = getDistance(lp1.get(i), lp2.get(i));
				sum += 1./(1 + d * d / d0sq);
			}
			value = (float) (sum/Lmin);
		}
		
		return new Tuple2<String, Float>(key.toString(), (float) value);
	}
	
	/**
	 * Get alignment for the two sequence. Object class casting.
	 * @param v1
	 * @param v2
	 * @param o
	 * @param e
	 * @return
	 */
	@SuppressWarnings("unchecked")
	private <T,K> Alignment<T> getAlignment(SequenceFeatureInterface<T> v1,SequenceFeatureInterface<K> v2,double o, double e) {
		return SmithWatermanGotoh.align(v1, (SequenceFeatureInterface<T>)v2, o, e);
	}
	
	/**
	 * Get the distance between two atoms
	 * @param a
	 * @param b
	 * @return
	 */
	public static double getDistance(Point3d a, Point3d b) {
		double x = a.x - b.x;
		double y = a.y - b.y;
		double z = a.z - b.z;
		double s  = x * x  + y * y + z * z;
		return Math.sqrt(s);
	}
	
	/**
	 * Rotate and transform based on the last alignment, in order to overlap two chains more accurately.
	 * @param SWAlignment	last alignment
	 * @param c1	chain1
	 * @param c2	chain2
	 * @param name2
	 * @param name
	 * @param it
	 * @return
	 * @throws FileNotFoundException
	 */
	private Alignment<?> rotateAlign(Alignment<?> SWAlignment, Point3d[] c1, Point3d[] c2, String name2, String name, int it) throws FileNotFoundException {
		// if last alignment fail, just return null
		if (SWAlignment != null && SWAlignment.getSequence1().length > 1) {
			Integer[] v1Order = SWAlignment.getSequence1();
			Integer[] v2Order = SWAlignment.getSequence2();
			ArrayList<Point3d> lp1 = new ArrayList<Point3d>();
			ArrayList<Point3d> lp2 = new ArrayList<Point3d>();
			for (int i = 0; i < v1Order.length; i++) {
				if (v1Order[i] != null && v2Order[i] != null) {
					lp1.add(c1[v1Order[i]]);
					lp2.add(c2[v2Order[i]]);
				}
			}
			
			Point3d[] p1 = new Point3d[lp1.size()];
			Point3d[] p2 = new Point3d[lp2.size()];
			for (int i = 0; i < p1.length; i++) {
				p1[i] = lp1.get(i);
				p2[i] = lp2.get(i);
			}
			// use qcp for rotate and transform
			SuperPositionQCP qcp = new SuperPositionQCP();
			qcp.set(p1, p2);
			Matrix4d m = qcp.getTransformationMatrix();
			SuperPositionQCP.transform(m, c2);
			
//			if (name.equals("1B0B.A") && name2.equals("4HHB.A")) {
//				PrintWriter writer = new PrintWriter("/Users/Chris/Documents/RCSB/Data/pdb/1B0B" + it + ".pdb");
//				Atom[] atoms = getCAAtoms(c2);
//				for (Atom a: atoms) {
//					writer.print(toPdb(a).replaceAll("nullGLU", " GLU"));
//				}
//				writer.close();
//			}
			
			SequenceFeatureInterface<Point3d> s1 = new Point3dFeature(c1);
			SequenceFeatureInterface<Point3d> s2 = new Point3dFeature(c2);
			return getAlignment(s1, s2, open, extend);
		}
		else 
			return null;
	}
	
	private static Atom[] getCAAtoms(Point3d[] points) {
		int gaps = 0;
		for (Point3d p: points) {
			if (p == null) {
				gaps++;
			}
		}
		Chain c = new ChainImpl();
		c.setChainID("A");		

		Atom[] atoms = new Atom[points.length-gaps];

		for (int i = 0, j = 0; i < points.length; i++) {
			if (points[i] != null) {
				atoms[j] = new AtomImpl();
				atoms[j].setName("CA");
				atoms[j].setPDBserial(i + 1);
				Group g = new AminoAcidImpl();
				g.setPDBName("GLU");
				g.addAtom(atoms[j]);
				g.setResidueNumber("A", i + 1, null);
				c.addGroup(g);

				atoms[j].setX(points[i].x);
				atoms[j].setY(points[i].y);
				atoms[j].setZ(points[i].z);
				j++;
			}
		}

		return atoms;
	}
	
	private static String toPdb(Atom a) {
		StringBuffer w = new StringBuffer();
		DecimalFormat d3 = (DecimalFormat)NumberFormat.getInstance(Locale.US);
		d3.setMaximumIntegerDigits(4);
		d3.setMinimumFractionDigits(3);
		d3.setMaximumFractionDigits(3);
		
		DecimalFormat d2 = (DecimalFormat)NumberFormat.getInstance(Locale.US);
		d2.setMaximumIntegerDigits(3);
		d2.setMinimumFractionDigits(2);
		d2.setMaximumFractionDigits(2);
		
		String newline = System.getProperty("line.separator");
		
		Group g = a.getGroup();

		GroupType type = g.getType() ;

		String record = "" ;
		if ( type.equals(GroupType.HETATM) ) {
			record = "HETATM";
		} else {
			record = "ATOM  ";
		}


		// format output ...
		//int groupsize  = g.size();
		String resName = g.getPDBName(); 
		String pdbcode = g.getResidueNumber().toString();
		//String line    = "" ;


		int    seri       = a.getPDBserial()        ;
		String serial     = alignRight(""+seri,5)   ;		
		String fullName   = formatAtomName(a);

		
		// System.out.println(" fullname: " + fullname + " : " + a.getAltLoc() + " : " + pdbcode);

		Character  altLoc = a.getAltLoc()           ;
		String resseq = "" ;
		if ( hasInsertionCode(pdbcode) )
			resseq     = alignRight(""+pdbcode,5);
		else
			resseq     = alignRight(""+pdbcode,4)+" ";
		String x          = alignRight(""+d3.format(a.getX()),8);
		String y          = alignRight(""+d3.format(a.getY()),8);
		String z          = alignRight(""+d3.format(a.getZ()),8);
		String occupancy  = alignRight(""+d2.format(a.getOccupancy()),6) ;
		String tempfactor = alignRight(""+d2.format(a.getTempFactor()),6);

		//System.out.println("fullname,size:" + fullname + " " + fullname.length());

		String leftResName = alignLeft(resName,3);

		StringBuffer s = new StringBuffer();
		s.append(record);
		s.append(serial);
		s.append(" ");
		s.append(fullName);
		s.append(altLoc);
		s.append(leftResName);
		s.append(" ");
		s.append(a.getGroup().getChain().getChainID());
		s.append(resseq);
		s.append("   ");
		s.append(x);
		s.append(y);
		s.append(z);
		s.append(occupancy);
		s.append(tempfactor);

		Element e = a.getElement();
		
		String eString = e.toString().toUpperCase();
		
		if ( e.equals(Element.R)) {
			eString = "    X";
		}
		w.append(String.format("%-76s%2s", s.toString(),eString));
		w.append(newline);
		
		
		return w.toString();
	}
	
	private static boolean hasInsertionCode(String pdbserial) {
		try {
			Integer.parseInt(pdbserial) ;
		} catch (NumberFormatException e) {
			return true ;
		}
		return false ;
	}
	
	private static String formatAtomName(Atom a) {
		String fullName = null;
		String name = a.getName();
		Element element = a.getElement();
		
		if (name.length()==4) 
			fullName = name;
		
		else if (name.length()==3) 
			fullName = " "+name;
				 
		else if (name.length()==2) {
			if (element == Element.C || element == Element.N || element == Element.O || element == Element.P || element == Element.S) 
				fullName = " "+name+" "; 
			else 
				fullName = name+"  ";
		}
		
		else if (name.length()==1) 
			fullName = " "+name+"  ";

		return fullName;
	}
	
	private static String alignRight(String input, int length){
		int n = input.length();
		if ( n >= length)
			return input;

		String spaces = "                           " ;
		int diff = length - n ;
		StringBuffer s = new StringBuffer();

		s.append(spaces.substring(0,diff));
		s.append(input);

		return s.toString();
	}
	
private static String alignLeft(String input, int length){
		if (input.length() >= length) {
			return input;
		}

		String spaces = "                           " ;
		input += spaces.substring(0, length - input.length() );
		return input;

	}
	
	@Override
	public void setCoords(Broadcast<List<Tuple2<String, Point3d[]>>> coords) {
		this.coords = coords;
	}
}
