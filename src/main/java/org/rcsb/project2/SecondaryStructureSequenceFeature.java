package org.rcsb.project2;

import scala.Tuple2;

public class SecondaryStructureSequenceFeature implements SequenceFeature<SecondaryStructureSequenceFeature, Byte> {

	byte[] feat;

	@SafeVarargs
	public SecondaryStructureSequenceFeature(int len, Tuple2<int[], int[]>... features) {
		feat = new byte[len];
		for (int k = 0; k < features.length; k++) {
			int[] start = features[k]._1;
			int[] end = features[k]._2;
			if (start.length != end.length)
				throw new IllegalArgumentException("Pair " + k + " array lengths not equal.");
			for (int i = 0; i < start.length; i++)
				for (int j = start[i]; j < end[i]; j++)
					feat[j] |= 1 << k;
		}
	}

	@Override
	public double similarity(SequenceFeature<SecondaryStructureSequenceFeature, Byte> sequence2, int i, int j) {
		int dif = (get(i) ^ sequence2.get(j));
		int val = 0;
		if (dif == 1 || dif == 2)
			val = 1;
		if (dif == 3)
			val = 2;
		return 1 / (1.0 + val);
	}

	@Override
	public boolean identity(SequenceFeature<SecondaryStructureSequenceFeature, Byte> sequence2, int i, int j) {
		return get(i) == sequence2.get(j);
	}

	@Override
	public Byte[] getSequence() {
		Byte[] o = new Byte[feat.length];
		for (int i = 0; i < feat.length; i++)
			o[i] = feat[i];
		return o;
	}

	@Override
	public Byte get(int index) {
		return feat[index];
	}

	@Override
	public int length() {
		return feat.length;
	}

	@Override
	public String toString(int index) {
		switch (get(index).byteValue()) {
		case 0:
			return "";
		case 1:
			return "A";
		case 2:
			return "B";
		case 3:
			return "AB";
		default:
			return "" + get(index);
		}
	}
}