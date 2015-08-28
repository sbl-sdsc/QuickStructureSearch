package org.rcsb.ProteinLigandInteractionSearch;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.RowFactory;
import scala.Tuple2;
/**
 * Create rows of features for the Java RDD
 * @author Hinna Shabir
 *
 */
public class HadoopToParqRow implements Function<Tuple2<String, Iterable<String>>, Row> {
	public Row call(Tuple2<String, Iterable<String>> tuple) throws Exception {
		ArrayList <String> v= new ArrayList <String>();
		String [] keys= tuple._1.split(",");
		for (String s:tuple._2) {
			v.add(s);
		}
		String[] pdbIds = v.toArray(new String[v.size()]);

		List<String> aminos=Arrays.asList("ARG", "HIS","LYS","ASP","GLU","SER","THR","ASN","GLN","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP");
		String index= new String();
		index= "Other";
		for(String amino:aminos){
			if (keys[6].equals(amino)){
				index=keys[6] ; // Create an index field that has the residue name of the protein i.e. 21 distinct values
				break;
			}
		}
		int size=keys.length;
		if (size!=13){
		System.out.println("Length: "+size);
		System.out.println(Arrays.toString(keys));
		System.out.println(Arrays.toString(pdbIds));
		}
		Integer dist= Integer.parseInt(keys[size-1]);
		return RowFactory.create(index,keys[0],keys[1],keys[2],keys[3],keys[4],keys[5],keys[6],keys[7],keys[8],keys[9],keys[10],keys[11],dist,pdbIds);
	}	
}