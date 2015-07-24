package org.rcsb.project2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.function.Predicate;

import javax.vecmath.Vector2d;

import scala.Tuple2;

public class SecondaryStructProjection {
	private final double X_DIFF = 5.0;
	private final double Y_DIFF = 5.0;
	private final int N;
	private Vector2d[] s, e;
	private Integer[] mapXS;
	private Integer[] mapYS;

	public SecondaryStructProjection(Vector2d[] s, Vector2d[] e) {
		if (s.length != e.length)
			throw new IllegalArgumentException("Lengths do not match!");
		this.s = s;
		this.e = e;
		this.N = s.length;
		this.mapXS = new Integer[N];
		this.mapYS = new Integer[N];
		for (int i = 0; i < N; i++) {
			mapXS[i] = i;
			mapYS[i] = i;
		}
		Arrays.sort(mapXS, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				double v = s[o1].x - s[o2].x;
				return v < 0 ? -1 : v == 0 ? 0 : 1;
			}
		});
		Arrays.sort(mapYS, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				double v = s[o1].y - s[o2].y;
				return v < 0 ? -1 : v == 0 ? 0 : 1;
			}
		});
		// System.out.println(Arrays.toString(s));
		// System.out.println(Arrays.toString(e));
		// System.out.println(Arrays.toString(mapXS));
		// System.out.println(Arrays.toString(mapYS));
	}

	public int length() {
		return N;
	}

	public Vector2d getStart(int i) {
		return s[i];
	}

	public Vector2d getEnd(int i) {
		return e[i];
	}

	public Tuple2<Vector2d, Vector2d> get(int i) {
		return new Tuple2<>(new Vector2d(s[i]), new Vector2d(e[i]));
	}

	public Tuple2<Integer, Double> getCloseTo(Tuple2<Vector2d, Vector2d> t) {
		return getCloseTo(t._1, t._2);
	}

	public Tuple2<Integer, Double> getCloseTo(Vector2d st, Vector2d en) {
		Tuple2<Integer, Double> NO_MATCH = new Tuple2<>(-1, -1.0);
		// System.out.println("finding " + st + " , " + en);
		List<Integer> candX = new ArrayList<>();
		int indX = search(N, new Predicate<Integer>() {
			@Override
			public boolean test(Integer t) {
				return st.x - X_DIFF > s[mapXS[t]].x;
			}
		});
		// System.out.println("indX: " + indX + ", " + mapXS[indX] + ", " + s[mapXS[indX]]);

		// do {
		// // System.out.println("x is " + s[mapXS[indX]].x);
		// candX.add(mapXS[indX]);
		// } while (++indX < N && s[mapXS[indX]].x < st.x + X_DIFF);

		while (indX < N && s[mapXS[indX]].x < st.x + X_DIFF) {
			candX.add(mapXS[indX++]);
		}
		// System.out.println("CANDX: " + candX);
		List<Integer> candY = new ArrayList<>();
		int indY = search(N, new Predicate<Integer>() {
			@Override
			public boolean test(Integer t) {
				return st.y - Y_DIFF > s[mapYS[t]].y;
			}
		});
		while (indY < N && s[mapYS[indY]].y < st.y + Y_DIFF) {
			candY.add(mapYS[indY++]);
		}
		// System.out.println("X: " + candX);
		// for (int i : candX) {
		// System.out.println("possible x : " + s[i].x);
		// }
		// System.out.println("Y: " + candY);
		// for (int i : candY) {
		// System.out.println("possible y : " + s[i].y);
		// }
		candX.retainAll(candY);
		// System.out.println("remaining candidates" + candX);
		// for (int i : candX) {
		// System.out.println("cand " + i + " : " + s[i]);
		// }
		if (candX.size() == 0)
			return NO_MATCH;
		int min = candX.get(0);
		double score = Double.MAX_VALUE;
		for (int i : candX) {
			double d = SecondaryStructTools.simil(st, en, s[i], e[i]);
			// System.out.println("Score for " + i + ", " + d);
			if (d < score) {
				score = d;
				min = i;
			}
		}
		// System.out.println(min + " Score: " + score);
		return score < 50 ? new Tuple2<>(min, score) : NO_MATCH;
	}

	private static int search(int len, Predicate<Integer> bigger) {
		int lo = 0;
		int hi = len - 1;
		int at = (lo + hi) / 2;
		while (hi - lo > 1) {
			if (bigger.test(at))
				lo = at;
			else
				hi = at;
			at = (lo + hi) / 2;
		}
		return bigger.test(lo) ? lo + 1 : lo;
	}
}
