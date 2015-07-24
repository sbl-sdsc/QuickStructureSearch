package org.rcsb.project2;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.StructureException;

public class SecondaryStructTest2 {

	private static final String NAME = "1NB5.J";

	public static void main(String[] args) {
		File f = new File("data/" + NAME + ".txt");
		SecondaryStruct s = null;
		if (f.exists())
			s = new SecondaryStruct(SecondaryStructTools.read(NAME));
		else {
			Point3d[] pts = null;
			try {
				pts = SecondaryStructTools.pull(NAME);
			}
			catch (IOException | StructureException e) {
				e.printStackTrace();
			}
			SecondaryStructTools.write(pts, NAME);
			s = new SecondaryStruct(pts);
		}
		SecondaryStructureSequenceFeature sf = s.getSequenceFeature();
		System.out.println("Start");
		System.out.println(NAME);
		try (Scanner scan = new Scanner(System.in)) {
			String in;
			while (!(in = scan.next()).equals("X")) {
				if (in.equals("g")) {
					int st = scan.nextInt();
					System.out.println(SecondaryStructTools.distsToString(s.getRange(st, scan.nextInt()), st));
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
					s.testPrint(s.getHelices());
				else if (in.equals("test1"))
					System.out.println("=(0,0,0)\t=" + s.normP + "*50\t=" + s.normX + "*50");
				else if (in.equals("test2"))
					SecondaryStruct.printProjection(s.getAlphaNormProjection((byte) 0b00000000));
			}
		}
		// sc.close();
	}
}
