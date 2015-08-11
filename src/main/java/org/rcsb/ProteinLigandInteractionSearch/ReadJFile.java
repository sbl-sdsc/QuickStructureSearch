package org.rcsb.ProteinLigandInteractionSearch;

import java.io.FileReader;
import java.util.Iterator;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

public class ReadJFile {

	public String ReadJ () {
		String res1=null;
		String res2= null;
		String atom1=null;
		String atom2=null;
		int distmin = 0;
		int distmax = 0;
		JSONParser parser = new JSONParser();
		try {
			Object obj = parser.parse(new FileReader("/Users/hina/query1.txt"));
			JSONObject jsonObject = (JSONObject) obj;
			JSONArray interactions = (JSONArray) jsonObject.get("interactions");
			Iterator<?> iterator = interactions.iterator();
			while (iterator.hasNext()) {
				JSONObject innerobj = (JSONObject) iterator.next();
//				System.out.println("res1 :"+innerobj.get("res1")+"res2 :"+innerobj.get("res2")+"distance:"+innerobj.get("distanceMin")+"-"+innerobj.get("distanceMax"));
				int site1 = Integer.parseInt((String) innerobj.get("res1"));
				int site2 = Integer.parseInt((String) innerobj.get("res2"));
				distmin =Integer.parseInt((String) innerobj.get("distanceMin"));
				distmax =Integer.parseInt((String) innerobj.get("distanceMax"));
				System.out.println(distmin);
				JSONArray sites= (JSONArray) jsonObject.get("sites");
				JSONObject innerObject1 =(JSONObject) sites.get(site1-1);
				JSONObject innerObject2 =(JSONObject) sites.get(site2-1);
				res1 =  (String) innerObject1.get("residueName");
				res2 =  (String) innerObject2.get("residueName");
				atom1 = (String) innerObject1.get("atomName");
				atom2 = (String) innerObject2.get("atomName");
				System.out.println(atom1);
				System.out.println(atom2);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		return (res1+"-"+res2+"-"+atom1+"-"+atom2+"-"+distmin+"-"+distmax);

	}
}