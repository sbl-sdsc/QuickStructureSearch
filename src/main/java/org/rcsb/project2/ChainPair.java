package org.rcsb.project2;


class ChainPair {
	private String n1, n2;
	private double tm;

	public String getN1() {
		return n1;
	}

	public void setN1(String n1) {
		this.n1 = n1;
	}

	public String getN2() {
		return n2;
	}

	public void setN2(String n2) {
		this.n2 = n2;
	}

	public double getTm() {
		return tm;
	}

	public void setTm(double tm) {
		this.tm = tm;
	}

	@Override
	public String toString() {
		return String.format("%s vs %s :: %.5f", getN1(), getN2(), getTm());
	}
}