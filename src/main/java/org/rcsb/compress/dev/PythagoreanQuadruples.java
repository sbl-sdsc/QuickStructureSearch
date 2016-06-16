package org.rcsb.compress.dev;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author Peter
 */
public class PythagoreanQuadruples {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
    	int dMin = 3600;
    	int dMax = 4000;
        int max = 100;
        int t = 0;

        long start = System.nanoTime();
        
        List<Integer[]> list = getQuadruples(dMin, dMax, max);
        
        long end = System.nanoTime();
        System.out.println("Time: " + (end - start) / 1E6 + " ms");
        System.out.println("Number of pythagorean triples t: " + t + " out of: " + max + " list.size: " + list.size());

        list = getReorderedQuadruples(dMin, dMax, max);
        for (Integer[] a : list) {
            System.out.println(Arrays.toString(a));
        }
    }

	public static List<Integer[]> getQuadruples(int dMin, int dMax, int max) {
		List<Integer[]> list = new ArrayList<Integer[]>();
  
        for (int m = 0; m < max; m++) {
            int mSq = m * m;
            for (int n = m; n < max; n++) {
                int nSq = n * n;
                int nmSq = mSq + nSq;
                if (nmSq > dMax) {
                    break;
                }
                for (int p = 1; p < max; p++) {
                    int d = p * p + nmSq;
                    if (d < dMin) {
                        continue;
                    }
                    if (d > dMax) {
                        break;
                    }
                    int a = 2 * m * p;
                    int b = 2 * n * p;
                    int c = p * p - nmSq;
                    if (c < 0) {
                        continue;
                    }

                    Integer[] v = new Integer[4];
                    v[0] = d;
                    v[1] = a;
                    v[2] = b;
                    v[3] = c;

                    Arrays.sort(v, Comparator.reverseOrder());
                    list.add(v);
                }
            }
        }
        Collections.sort(list, QuadrupleComparator);
		return list;
	}
    
	public static List<Integer[]> getReorderedQuadruples(int dMin, int dMax, int max) {
		List<Integer[]> list = getQuadruples(dMin, dMax, max);
		List<Integer[]> reordered = new ArrayList<Integer[]>(list.size());
		if (list.size() % 2 == 0) {
			int middle = list.size()/2;
			for (int i = 0; i < middle; i++) {
				reordered.add(list.get(middle-i-1));
				reordered.add(list.get(middle+i));
			}
		}
		return reordered;
	}
	
    public static Comparator<Integer[]> QuadrupleComparator 
                          = new Comparator<Integer[]>() {

	    public int compare(Integer[] quadruple1, Integer[] quadruple2) {
	      
	      //ascending order
              int delta = quadruple1[0] - quadruple2[0];
              if (delta != 0) {
                  return delta;
              }
              
	      return quadruple1[1] - quadruple2[1];
	    }

	};
}