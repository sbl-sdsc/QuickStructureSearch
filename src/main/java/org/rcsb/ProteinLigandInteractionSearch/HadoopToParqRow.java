package org.rcsb.ProteinLigandInteractionSearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.spark.api.java.function.Function;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.RowFactory;

import scala.Tuple2;



public class HadoopToParqRow implements Function<Tuple2<String, Iterable<String>>, Row> {

	private static final long serialVersionUID = 1L;
	@Override
	public Row call(Tuple2<String, Iterable<String>> tuple) throws Exception {

		ArrayList <String> v= new ArrayList <String>();
		String [] keys= tuple._1.split("-");
		for (String s:tuple._2) {
			v.add(s);
		}
		
		String[] pdbIds = v.toArray(new String[v.size()]);
		// System.out.println("# of PdbIds: " + pdbIds.length);
		 
	/*	for (int i=0; i<pdbIds.length ;i++){
	       System.out.println(i);
	       System.out.println("pdbIds :"+ pdbIds[i]);
		}*/
		int flag=0;
		List<String> aminos=Arrays.asList("ARG", "HIS","LYS","ASP","GLU","SER","THR","ASN","GLN","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP");
		String index= new String();
		for(String amino:aminos){
			if (keys[0].equals(amino)){
				flag=1;
				break;
			}
		}
		if (flag==1){
			index=keys[0] ;
		}
		else{
			index= "Other";
		}
		return RowFactory.create(index,keys[0],keys[1],keys[2],keys[3],keys[4],keys[5],Integer.parseInt(keys[6]),pdbIds);

	}
	
}