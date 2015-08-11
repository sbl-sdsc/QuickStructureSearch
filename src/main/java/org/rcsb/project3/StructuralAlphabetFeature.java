package org.rcsb.project3;

import java.io.Serializable;

import javax.vecmath.Point3d;

/**
 * This class implements the SequenceFeatureInterface.
 * It is used for the calculation for StructuralAlphabetFingerprint
 * 
 * @author Chris Li
 */
public class StructuralAlphabetFeature implements SequenceFeatureInterface<String>, Serializable {

	private static final long serialVersionUID = 1L;
	private String[] AlphabetSequence;
	private Point3d[] coords;
	// Some setting for calculate similarity score 
	private double match = 1;
    private double mismatch = -1;
    private double gap = -1;
    // score matrix for checking similarity
    private static double[][] scoreMatrix = {
    		{ 516, -59, 113, -105, -411,-177, -27,-361,  47,-103,-644,-259, -599, -372,-124, -83 },
    		{ -59, 541,-146, -210, -155,-310, -97,  90, 182,-128, -30,  29, -745, -242,-165,  22 },
    		{ 113,-146, 360,  -14, -333,-240,  49,-438,-269,-282,-688,-682, -608, -455,-147,   6 },  
    		{-105,-210, -14,  221,    5,-131,-349,-278,-253,-173,-585,-670,-1573,-1048,-691,-497 },
    		{-411,-155,-333,    5,  520, 185, 186, 138,-378, -70,-112,-514,-1136, -469,-617,-632 },
    		{-177,-310,-240, -131,  185, 459, -99, -45,-445,  83,-214, -88, -547, -629,-406,-552 },
    		{ -27, -97,  49, -349,  186, -99, 665, -99, -89,-118,-409,-138, -124,  172, 128, 254 },
    		{-361,  90,-438, -278,  138, -45, -99, 632,-205, 316, 192,-108, -712, -359,  95,-399 },
    		{  47, 182,-269, -253, -378,-445, -89,-205, 696, 186,   8,  15, -709, -269,-169, 226 },
    		{-103,-128,-282, -173,  -70,  83,-118, 316, 186, 768, 196,   5, -398, -340,-117,-104 },
    		{-644, -30,-688, -585, -112,-214,-409, 192,   8, 196, 568, -65, -270, -231,-471,-382 },
    		{-259,  29,-682, -670, -514, -88,-138,-108,  15,   5, -65, 533, -131,    8, -11,-316 },
    		{-599,-745,-608,-1573,-1136,-547,-124,-712,-709,-398,-270,-131,  241,   -4,-190,-155 },
    		{-372,-242,-455,-1048, -469,-629, 172,-359,-269,-340,-231,   8,   -4,  703,  88, 146 },
    		{-124,-165,-147, -691, -617,-406, 128,  95,-169,-117,-471, -11, -190,   88, 716,  58 },
    		{ -83,  22,   6, -497, -632,-552, 254,-399, 226,-104,-382,-316, -155,  146,  58, 609 },
		};
    private static String blockName = "abcdefghijklmnop";
	
    /**
     * Constructor that will store a String array
     * @param AngleSequence
     */
	public StructuralAlphabetFeature(String[] AlphabetSequence) {
		this.AlphabetSequence = AlphabetSequence;
	}
	
	public StructuralAlphabetFeature(String[] AlphabetSequence, Point3d[] coords) {
		this.AlphabetSequence = AlphabetSequence;
		this.coords = coords;
	}
	
	/**
	 * Constructor that will store a double array of angle and update the settings
	 * @param AngleSequence
	 * @param diff
	 * @param gap
	 * @param match
	 * @param mismatch
	 */
	public StructuralAlphabetFeature(String[] AlphabetSequence, double gap, double match, double mismatch) {
		this.AlphabetSequence = AlphabetSequence;
		this.gap = gap;
		this.match = match;
		this.mismatch = mismatch;
	}
		
	@Override
	public double similarity(SequenceFeatureInterface<String> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return gap;
		}
		// check similarity
		else if (this.get(i).equals(sequence2.get(j)))
			return match;
		else {
			int index1 = blockName.indexOf(this.get(i));
			int index2 = blockName.indexOf(sequence2.get(j));
			if (index1 < 0 || index2 < 0)
				return mismatch;
			double score = (scoreMatrix[index1][index2] + 1573) / (768 + 1573);
			return score;
		}
	}

	@Override
	public boolean identity(SequenceFeatureInterface<String> sequence2, int i, int j) {
		// check NaN as gap
		if (this.get(i) == null || sequence2.get(j) == null){
			return false;
		}
		// check identity
		else if (this.get(i).equals(sequence2.get(j)))
			return true;
		else 
			return false;
	}

	@Override
	public String get(int index) {
		return AlphabetSequence[index];
	}

	@Override
	public int length() {
		return AlphabetSequence.length;
	}

	@Override
	public String[] getSequence() {
		return AlphabetSequence;
	}

	@Override
	public String toString(int index) {
		return this.get(index);
	}

	@Override
	public double todouble(int index) {
		return 0;
	}

	@Override
	public Point3d[] getCoords() {
		return coords;
	}
}
