package org.rcsb.project10;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.rcsb.mmtf.api.StructureDataInterface;

import scala.Tuple2;

/**
 * A mapper from {@link StructureDataInterface} to the {@link Point3d}[] of the calpha coordinates.
 * @author Anthony Bradley
 * @author Peter W. Rose
 *
 */
public class GappedSegmentGenerator implements PairFlatMapFunction<Tuple2<String,StructureDataInterface>,String, WritableSegment> {
	private static final long serialVersionUID = 1535046883440398492L;
	
	/** Minimum number of C alpha atoms with 3D coordinates */
	private int minSequenceLength;

	/**
	 * 
	 * @param fragmentLength the length of each fragment. Null means take each Chain as 
	 * as single fragment.
	 */
	public GappedSegmentGenerator(int minSequenceLength) {
		this.minSequenceLength = minSequenceLength;
	}

	@Override
	public Iterable<Tuple2<String, WritableSegment>> call(Tuple2<String, StructureDataInterface> t) throws Exception {
		StructureDataInterface structureDataInterface = t._2;
		
		List<Tuple2<String, WritableSegment>> outList = new ArrayList<>();
		
		String pdbId = structureDataInterface.getStructureId();
		
		// extract chain information
		boolean[] polymer = isPolymer(structureDataInterface);
		String[] sequences = getSequences(structureDataInterface);
		int[] groupSequenceIndices = structureDataInterface.getGroupSequenceIndices();
		
		int atomCounter = 0;
		int groupCounter = 0;

		// loop through chains
		for (int i = 0; i < structureDataInterface.getNumChains(); i++) {
			String chainId = pdbId + "." + structureDataInterface.getChainIds()[i];

			Point3d[] coords = null;
			String sequence = null;

			if (polymer[i]){			
				sequence = sequences[i];
				if (sequence.length() < minSequenceLength) {
					polymer[i] = false;
				}
				
				coords = new Point3d[sequence.length()];
				for (int j = 0; j < coords.length; j++) {
					coords[j] = null;
				}
			}

			int calphaCounter = 0;

			for (int groupId=0; groupId < structureDataInterface.getGroupsPerChain()[i]; groupId++){
				int groupType = structureDataInterface.getGroupTypeIndices()[groupCounter];

				// save C alpha coordinates
				if (polymer[i]) {
					Point3d point3d = getCalpha(structureDataInterface, groupType, atomCounter);
					// groupSequenceIndices can be -1, for example for heterogeneity: 4KL8.B
					if (point3d != null && groupSequenceIndices[groupCounter] >= 0) {
						coords[groupSequenceIndices[groupCounter]] = point3d;
						calphaCounter++;
					}
				}
				atomCounter += structureDataInterface.getNumAtomsInGroup(groupType);
				groupCounter++;	
			}

			if (calphaCounter >= minSequenceLength) {
				outList.add(new Tuple2<String, WritableSegment>(chainId, new WritableSegment(sequence, coords)));
			}
		}
		return outList;
	}
	
	/**
	 * Returns a boolean array that indicated which chains are polymers.
	 * @param structureDataInterface the input {@link StructureDataInterface}
	 * @return a boolean array, if true, chain is a polymer
	 */
	private boolean[] isPolymer(StructureDataInterface structureDataInterface) {
		boolean[] polymer = new boolean[structureDataInterface.getNumChains()];

		for(int i = 0; i < structureDataInterface.getNumEntities(); i++) {
			String type = structureDataInterface.getEntityType(i);
			if (type.equals("polymer")) {
				for (int j : structureDataInterface.getEntityChainIndexList(i)) {
					polymer[j] = true;
				}
			}
		}
		return polymer;
	}

	/**
	 * Gets the polymer one-letter sequence for each chain
	 * @param structureDataInterface the input {@link StructureDataInterface}
	 * @return array of String containing the polymer sequence for each chain
	 */
	private String[] getSequences(StructureDataInterface structureDataInterface) {
		String[] sequences = new String[structureDataInterface.getNumChains()];
		
		for (int i = 0; i < structureDataInterface.getNumEntities(); i++) {
			String sequence = structureDataInterface.getEntitySequence(i);
			for (int j : structureDataInterface.getEntityChainIndexList(i)) {
				sequences[j] = sequence;
			}
		}
		return sequences;
	}


	/**
	 * Gets the calpha coordinates of a group as a {@link Point3d}.
	 * @param structureDataInterface the {@link StructureDataInterface} to read
	 * @param groupType the integer specifying the grouptype
	 * @param atomCounter the atom count at the start of this group
	 * @return the point3d object specifying the calpha of this point
	 */
	private Point3d getCalpha(StructureDataInterface structureDataInterface, int groupType, int atomCounter) {
		for (int i = 0; i < structureDataInterface.getNumAtomsInGroup(groupType); i++){
			if (structureDataInterface.getGroupAtomNames(groupType)[i].equals("CA")){
				Point3d point3d = new Point3d();
				point3d.x = structureDataInterface.getxCoords()[atomCounter+i];
				point3d.y = structureDataInterface.getyCoords()[atomCounter+i]; 
				point3d.z = structureDataInterface.getzCoords()[atomCounter+i];
				return point3d;
			}
		}
		return null;
	}
}
