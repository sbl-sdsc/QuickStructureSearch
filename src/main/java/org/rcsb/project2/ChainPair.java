package org.rcsb.project2;

/**
 * Class for a pair of proteins that needs to be compared
 * 
 * @author Kevin Wu
 *
 */
class ChainPair {
	private String name1, name2;
	private double tmScore;

	/**
	 * 
	 * @return String of the first protein chain
	 */
	public String getN1() {
		return name1;
	}

	/**
	 * Sets the first protein chain name
	 * 
	 * @param n1
	 *            name of first protein chain
	 */
	public void setN1(String n1) {
		this.name1 = n1;
	}

	/**
	 * 
	 * @return String of the second protein chain
	 */
	public String getN2() {
		return name2;
	}

	/**
	 * Sets the second protein chain name
	 * 
	 * @param n2
	 *            name of second protein chain
	 */
	public void setN2(String n2) {
		this.name2 = n2;
	}

	/**
	 * 
	 * @return TM Score of the alignment of the two protein chains
	 */
	public double getTm() {
		return tmScore;
	}

	/**
	 * Sets the TM Score for the alignment of the two protein chains
	 * 
	 * @param tm
	 */
	public void setTm(double tm) {
		this.tmScore = tm;
	}

	/**
	 * toString method for printing
	 */
	@Override
	public String toString() {
		return String.format("%s vs %s :: %.5f", getN1(), getN2(), getTm());
	}
}