package org.rcsb.project2;

import java.util.Arrays;

import scala.Tuple2;

/**
 * For merging alpha helices and beta strands
 * 
 * @author Kevin Wu
 *
 */
public class MergeStruct extends SecondaryStructFeature {

	public MergeStruct(SecondaryStruct s) {
		this.pts = s.getPoints();
		this.normP = s.normP;
		this.normX = s.normX;
		this.normC = s.normC;
		int[] start, end;
		int N = s.getAlphaLength() + s.getBetaLength();
		Tuple2<int[], int[]> ta = s.getAlpha().getFeatures();
		Tuple2<int[], int[]> tb = s.getBeta().getFeatures();
		int la = s.getAlphaLength();
		start = Arrays.copyOf(ta._1, N);
		end = Arrays.copyOf(ta._2, N);
		for (int i = la; i < N; i++) {
			start[i] = tb._1[i - la];
			end[i] = tb._2[i - la];
		}
		features = new Tuple2<>(start, end);
		projections = new SecondaryStructProjection[N];
	}
}
