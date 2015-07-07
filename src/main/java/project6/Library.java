package project6;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.hadoop.io.ArrayWritable;
import org.apache.hadoop.io.Text;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.rcsb.hadoop.io.Demo;
import org.rcsb.hadoop.io.HadoopToSimpleChainMapper;
import org.rcsb.hadoop.io.SimplePolymerChain;
import org.rcsb.structuralAlignment.SuperPositionQCP;
import org.rcsb.structuralSimilarity.GapFilter;
import org.rcsb.structuralSimilarity.LengthFilter;

import scala.Tuple2;

/**
 * 
 * This class creates a library of unique fragments categorized by length.
 * The only data to be output is a key for each fragment and an arraywritable
 * of length ranges and x0,y0,z0,x1,y1,z1,... coordinates.\
 * For now, though, I'm getting it to output length as a string and the points
 * as a point3d[].
 * 
 * @author Grant Summers
 *
 */

public class Library
{
	private static int NUM_THREADS = 4;
	private static int NUM_TASKS_PER_THREAD = 3;
	
	// make the data structure to store unique fragments
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException
	{
		 // arguments
		 String path = "/Users/grantsummers/Desktop/School/Internship/main/QuickStructureSearch/src/main/java/org/rcsb/project1/protein_chains_All_20150629_002251.seq";
//		 String outputfilename = "output.txt";
		 int length = 8;
		 
		// Spark setup
		SparkConf conf = new SparkConf().setMaster("local[" + NUM_THREADS + "]")
				.setAppName(Demo.class.getSimpleName())
				.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
		JavaSparkContext sc = new JavaSparkContext(conf);
		
		// start time
//		long start = System.nanoTime();
		
		// map sequence file to pairs (id, points)
		List<Tuple2<String, Point3d[]>> chains = sc  // CHANGE THIS TO SIMPLEPOLYMERCHAIN
				.sequenceFile(path, Text.class, ArrayWritable.class,NUM_THREADS*NUM_TASKS_PER_THREAD)
				.sample(false, 0.0003, 123)
				.mapToPair(new HadoopToSimpleChainMapper()) // convert input to <pdbId.chainId, SimplePolymerChain> pairs
				.filter(t -> t._2.isProtein())
				.mapToPair(t -> new Tuple2<String, Point3d[]>(t._1, t._2.getCoordinates()))
				.filter(new GapFilter(0, 0)) // keep protein chains with gap size <= 0 and 0 gaps
				.filter(new LengthFilter(50,500)) // keep protein chains with 50 - 500 residues
				.collect();		
		
		List<Tuple2<String, Point3d[]>> lib = new ArrayList<>();
		
		// index for where to place the fragments in lib
		int index = 0;
		
		// boolean for whether or not to add a fragment to the library
		boolean bool = true;
		
		// instantiate the SPQCP object
		SuperPositionQCP qcp = new SuperPositionQCP(true);
		
		for (Tuple2<String, Point3d[]> t: chains) {
			for(int star=0; star<t._2.length-length; star++)
			{
//				System.out.println(t._1 + "." + star + ": " + Arrays.toString(Arrays.copyOfRange(t._2, star, star + length)));

				// Create a Tuple2 for each fragment
				Tuple2<String, Point3d[]> tup = new Tuple2<String, Point3d[]>(t._1 + "." + star, Arrays.copyOfRange(t._2, star, star+length));
				
				// center each fragment								  //
				SuperPositionQCP.center(tup._2);					 //
																	//
																   //
																  //
				if(!lib.isEmpty()){								 //
					for(Tuple2<String, Point3d[]> l: lib){		//
						if(l._2 != null && tup._2 != null){	   //
							qcp.set(l._2, tup._2);			  //
							if(qcp != null){				 //
								double q = qcp.getRmsd();	// This line gives a null pointer exception
								if(q<1){				  //\\
									bool = false;		   //\\
								}							//\\
							}								 //\\
						}									  //\\
						else{								   //\\
							bool = false;						//\\
						}										 //\\
					}											  //\\
				}												   //\\
				if(bool == true){									//\\
					tup.copy(index, tup._2);
					lib.add(index, tup);
					System.out.println(index + ": " + Arrays.toString(lib.get(index)._2));
					index++;
				}
				bool = true;
			}
		}
		sc.close();
		
		// prints all fragments in lib
//		for(Tuple2<String, Point3d[]> l: lib){
//			System.out.println(l._1 + Arrays.toString(l._2));
//		}
		
		// Write the lib list to a text or csv file
		PrintWriter writer = new PrintWriter("library.txt", "UTF-8");
		for(Tuple2<String, Point3d[]> l: lib){
			writer.println(l._1 + Arrays.toString(l._2));
		}
		writer.close();
		
		// prints time
//		System.out.println("Time: " + (System.nanoTime() - start)/1E9 + " sec.");
	}
	
	public String lengthy(Point3d[] p){
		double i = Math.abs(p[p.length-1].distance(p[0]));
		i /= 2;
		int base = (int) i;
		base *= 2;
		int top = base + 2;
		return Integer.toString(base) + " - " + Integer.toString(top);
	}
}
