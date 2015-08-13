package org.rcsb.project2;

/**
 * Interface for comparing secondary structures
 * 
 * @author Kevin Wu
 *
 */
public interface StructureCompare {
	public SecondaryStructProjection getNormProjection(byte b);

	public SecondaryStructProjection getProjection(int i);

	public int length();
}
