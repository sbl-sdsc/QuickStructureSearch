package org.rcsb.REHS;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;

import javax.vecmath.Point3d;

import org.apache.commons.lang.StringUtils;
import org.apache.hadoop.io.Text;
import org.apache.spark.Accumulable;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.project10.WritableSegment;

/**
 * This class generates a library of unique fragments from a protein chain file
 * when given a fragment size and and rmsd threshold (used to see whether 
 * or not a fragment is unique compared to the others). The library is stored as a text file.
 * This class also contains a static method to read the library into a list, 
 * which the LibraryFingerprint class can use.
 * 
 * @author Emilia Copic
 * @author Varkey Alumootil
 * @author Peter Rose
 *
 */
public class ArchLibGenerator {
	public static void main(String[] args) throws IOException {
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		System.out.println("ArchLibGen");

		String chainFile = args[0];
		System.out.println("Chain file    : " + chainFile);

		long start = System.nanoTime();

		int fragmentSize = 8;
		double rmsdThreshold = 2.05;

		getAllFragments(chainFile, fragmentSize, rmsdThreshold, args[1] + "_" + timeStamp);

		long end = System.nanoTime();

		System.out.println("Time          : " + (end - start) / 1E9 + " seconds");
	}
/**
 * Generates library, stores it in text file
 * @param chainFile file with all of the proteins
 * @param fragmentSize number of 3D points each fragment will be
 * @param rmsdThreshold when a fragment is compared to one in the library and the rmsd is below this threshold, it will not be considered unique and therefore not be added to the library
 * @param filePath library text file will be stored here
 * @throws IOException
 */
	private static void getAllFragments(String chainFile, int fragmentSize, double rmsdThreshold, String filePath)
			throws IOException {

		// setup spark
		SparkConf conf = new SparkConf().setMaster("local[*]").setAppName("ArchLibGenerator").set("spark.serializer",
				"org.apache.spark.serializer.KryoSerializer");

		JavaSparkContext sc = new JavaSparkContext(conf);

		// read protein chains and cut into fragments (sliding window approach)
		JavaRDD<Point3d[]> fragments = sc.sequenceFile(chainFile, Text.class, WritableSegment.class)
				.map(t -> t._2.getCoordinates()) // get the coordinates of the
													// protein chains
				.repartition(1) // create a single partition to generate a
								// single fragment library (this cannot be done
								// in parallel!)
				.flatMap(new FlatMapToFragments(fragmentSize)); // flatmap to
																// fragments

		// Create library of fragment archetypes
		List<Point3d[]> prelib = new ArrayList<>();
		List<Point3d[]> lib = Collections.synchronizedList(prelib);
		Accumulable<List<Point3d[]>, Point3d[]> accLib = new Accumulable<>(lib, new AccumuableList(rmsdThreshold));

		fragments.foreach(t -> accLib.add(t));

		List<Point3d[]> archetypes = accLib.value();

		System.out.println(archetypes.size());

		// save to text file for later use

		FileWriter writer = new FileWriter(filePath);
		for (Point3d[] fragment : archetypes) {
			String frag = Arrays.toString(fragment);
			writer.write(frag);

		}
		writer.close();

		sc.stop();
		sc.close();
	}
/**
 * Reads text file genreated by getAllFragments and returns the library as a list
 * @param filePath gets text file from this location
 * @return libray as a List
 * @throws FileNotFoundException
 */
	public static List<Point3d[]> readLibraryFromFile(String filePath) throws FileNotFoundException {
		// thank you for nO HELP VARKEY
		Scanner s = new Scanner(new File(filePath));
		String all = "";
		ArrayList<Point3d[]> lib = new ArrayList<>();
		while (s.hasNext()) {
			all += s.next();
		}
		// split into list
		String[] a = all.split("]");
		String[] b = new String[a.length];
		// get rid of "]"
		int fragLength = StringUtils.countMatches(a[0], ")");
		for (int i = 0; i < b.length; i++) {
			b[i] = a[i].replace("[", "");
			b[i] = b[i].replace("(", "");
			b[i] = b[i].replace(")", "");
			b[i] = b[i].replace(",", " ");
		}

		for (int i = 0; i < b.length; i++) {
			String pointArray = b[i];
			Scanner scanString = new Scanner(pointArray);
			Point3d[] frag = new Point3d[fragLength];
			for (int j = 0; j < fragLength; j++) {
				Point3d point3d = new Point3d(scanString.nextDouble(), scanString.nextDouble(),
						scanString.nextDouble());
				frag[j] = point3d;
			}
			scanString.close();
			lib.add(frag);
		}

		s.close();
		return lib;
	}

}
