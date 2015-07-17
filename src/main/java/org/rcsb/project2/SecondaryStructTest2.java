package org.rcsb.project2;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.StructureException;

public class SecondaryStructTest2 {

	private static final String NAME = "4FMW.B";

	public static void main(String[] args) {
		File f = new File(NAME + ".txt");
		SecondaryStruct s = null;
		if (f.exists())
			s = new SecondaryStruct(SecondaryStruct.read(NAME));
		else {
			Point3d[] pts = null;
			try {
				pts = SecondaryStruct.pull(NAME);
			}
			catch (IOException | StructureException e) {
				e.printStackTrace();
			}
			SecondaryStruct.write(pts, NAME);
			s = new SecondaryStruct(pts);
		}
		SecondaryStructureSequenceFeature sf = s.getSequenceFeature();
		System.out.println("Start");
		try (Scanner scan = new Scanner(System.in)) {
			String in;
			while (!(in = scan.next()).equals("X")) {
				if (in.equals("g")) {
					int st = scan.nextInt();
					System.out.println(SecondaryStruct.distsToString(s.getRange(st, scan.nextInt()), st));
				}
				else if (in.equals("a"))
					s.printHelices();
				else if (in.equals("b"))
					s.printStrands();
				else if (in.equals("l"))
					System.out.println(s.length());
				else if (in.equals("c"))
					s.printPoints();
				else if (in.equals("sf"))
					for (int i = 0; i < s.length(); i++)
						System.out.println((i + 1) + ":\t" + sf.toString(i));
				else if (in.equals("test"))
					s.printBetaProjection(0);
			}
		}
		// sc.close();
	}
}
