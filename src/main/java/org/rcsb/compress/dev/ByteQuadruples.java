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
public class ByteQuadruples {
	private static int DISTANCE_MIN = 3672;
	private static int DISTANCE_MAX = 3928;
//	private static int DISTANCE_MIN = 3600;
//	private static int DISTANCE_MAX = 4000;
//	private static int BIN_SIZE = 256;
	private static int BIN_SIZE = 128;

	private List<Integer[]> quadruples;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
       List<Integer[]> quadruples = getNewQuadruples(DISTANCE_MIN, DISTANCE_MAX);
       System.out.println("quadruples: " + quadruples.size());
       for (Integer[] q: quadruples) {
    	   System.out.println(Arrays.toString(q) + ": " + q[0]*BIN_SIZE + "," + q[1]*BIN_SIZE + "," + q[2]*BIN_SIZE);
       }
    }
    
    public ByteQuadruples(int minDistance, int maxDistance) {
    	quadruples = getNewQuadruples(minDistance, maxDistance);
    }
    
	public int[] applyTransformation(int[] coords) {
		int index = -1;
		int dSq = Integer.MAX_VALUE;
		int dxmin = coords[0];
		int dymin = coords[1];
		int dzmin = coords[2];

		for (int i = 0; i < quadruples.size(); i++) {
			Integer[] quad = quadruples.get(i);

			int dx = quad[0]*BIN_SIZE - coords[0];
//			if (dx < 0) continue;
			int dy = quad[1]*BIN_SIZE - coords[1];
//			if (dy < 0) continue;
			int dz = quad[2]*BIN_SIZE - coords[2];
//			if (dz < 0) continue;

//			int dTemp = dx*dx + dy*dy + dz*dz;
			int dTemp = Math.abs(dx)+ Math.abs(dy) + Math.abs(dz);

			if (dTemp < dSq) {
				dSq = dTemp;
				index = i;
				dxmin = dx;
				dymin = dy;
				dzmin = dz;
			}
		}
		
		int[] results = {index, dxmin, dymin, dzmin};
		return results;
	}
	
    public List<Integer[]> getQuadruples() {
    	return this.quadruples;
    }
    
	private static List<Integer[]> getNewQuadruples(int dMin, int dMax) {
		List<Integer[]> list = new ArrayList<Integer[]>();
		int dMinSq = dMin*dMin;
		int dMaxSq = dMax*dMax;
		int median = (dMin+dMax)/2;
		int dMedSq = median * median;
		System.out.println("median: " + median);
		
		int iMax = (dMax / BIN_SIZE);
		if (iMax * BIN_SIZE < dMax) {
			iMax++;
		}
		
		int bSq = BIN_SIZE * BIN_SIZE;

        for (int i = -iMax; i <= iMax; i++) {
            int iSq = i*i*bSq;
            for (int j = -iMax; j <= iMax; j++) {
                int jSq = j*j*bSq;
                for (int k = -iMax; k <= iMax; k++) {
                    int dSq = iSq + jSq + k * k * bSq;
  //                  System.out.println("dSq: " + dSq);
                    if (dSq < dMinSq) {
                        continue;
                    }
                    if (dSq > dMaxSq) {
                        continue;
                    }
                  
                    Integer[] v = new Integer[4];
                    v[0] = i;
                    v[1] = j;
                    v[2] = k;
 //                   v[3] = Math.abs(dMedSq - dSq);
                    v[3] = (int) Math.abs(median - Math.sqrt(dSq));

  //                  Arrays.sort(v, Comparator.reverseOrder());
                    list.add(v);
                }
            }
        }
 //       Collections.sort(list, QuadrupleComparator);
//        for (Integer[] l: list) {
//        	System.out.println(Arrays.toString(l));
//        }
		return list;
	}
	
    private static Comparator<Integer[]> QuadrupleComparator 
                          = new Comparator<Integer[]>() {

	    public int compare(Integer[] quadruple1, Integer[] quadruple2) {
	      
	      //ascending order
              return quadruple1[3] - quadruple2[3];
	    }

	};
	
	public static int getIndex(int[] coordinates) {
		int[] unsignedCoords = new int[3];
		unsignedCoords[0] = Math.abs(coordinates[0]);
		unsignedCoords[1] = Math.abs(coordinates[1]);
		unsignedCoords[2] = Math.abs(coordinates[2]);
		
		int min = coordMin(unsignedCoords);
		int max = coordMax(unsignedCoords);

		byte index = 0;
		if (unsignedCoords[0] > min && unsignedCoords[0] < max) {
			index |= 1;
		}
		if (unsignedCoords[0] == max) {
			index |= 1 << 1;
		}
		if (unsignedCoords[1] > unsignedCoords[2]) {
			index |= 1 << 2;
		}
		if (coordinates[0] < 0) {
			index |= 1 << 3;
		} 
		if (coordinates[1] < 0) {
			index |= 1 << 4;
		}
		if (coordinates[2] < 0) {
			index |= 1 << 5;
		}
		
		return index;
	}
	
	private static int coordMin(int[] values) {
		return Math.min(values[0],  Math.min(values[1], values[2]));
	}
	
	private static int coordMax(int[] values) {
		return Math.max(values[0],  Math.max(values[1], values[2]));
	}
}