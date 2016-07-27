package org.rcsb.projectec;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import javax.vecmath.Point3d;

public class ReadLibTester {
	public static void main(String[] args) throws FileNotFoundException {
	List<Point3d[]> lib = ArchLibGeneratorPR.readLibraryFromFile(args[0]);
//		for (int i = 0; i<10 ; i++) {
//			System.out.println(Arrays.toString(lib.get(i)));
//		}
	ArchLibGeneratorPR b = new ArchLibGeneratorPR();
		double[][] a = b.getRmsdArray(lib);
		System.out.println(Arrays.toString(a[0]));
		System.out.println(Arrays.toString(a[1]));
		System.out.println(Arrays.toString(a[2]));
	}
}
