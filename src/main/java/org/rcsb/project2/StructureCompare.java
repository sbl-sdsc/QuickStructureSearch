package org.rcsb.project2;

public interface StructureCompare {
	public SecondaryStructProjection getNormProjection(byte b);

	public SecondaryStructProjection getProjection(int i);

	public int length();
}
