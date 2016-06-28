package org.rcsb.projectva;
import java.util.Map;
import java.util.Map.Entry;
import java.io.ByteArrayOutputStream;
import java.util.zip.GZIPOutputStream;


public class NormalizedCompressionDistance {
	
//	public static<T> float distance(Map<T, Integer> featureCounts1, Map<T,Integer> featureCounts2)
//	{
//		int sum1 = 0;
//		for (Entry<T, Integer> entry: featureCounts1.entrySet()) {
//			sum1 += entry.getValue();
//		}
//		
//		int sum2 = 0;
//		for (Entry<T, Integer> entry: featureCounts1.entrySet()) {
//			sum2 += entry.getValue();
//		}
//		
//		int[] s = new int[sum1];
//		int i = 0;
//		
//		for (Entry<T, Integer> entry: featureCounts1.entrySet()) {
//			for (int j = 0; j < entry.getValue(); j++) {
//				s[i++] = featureCounts1.get(entry.getKey());
//			}
//		}
//		
//		int[] t = new int[sum2];
//		i = 0;
//		
//		for (Entry<T, Integer> entry: featureCounts2.entrySet()) {
//			for (int j = 0; j < entry.getValue(); j++) {
//				t[i++] = featureCounts2.get(entry.getKey());
//			}
//		}
//		
//		return (float)distance(s,t);
//	}
	
	public static double distance(int[] s, int[] t)
	{
	    if (s.length == 0) return t.length;
	    if (t.length == 0) return s.length;

        int slen = length(s);
        int tlen = length(t);
        int clen = length(combine(s,t));
        int clen2 = length(combine(t,s));
        clen = Math.min(clen, clen2);
 //       System.out.println("s,t,c; " + slen + "," + tlen + "," + clen + ","  + clen2);
//        System.out.println((clen - Math.min(slen, tlen))/(double)Math.max(slen,  tlen));
//        System.out.println((clen2 - Math.min(slen, tlen))/(double)Math.max(slen,  tlen));
        
		return 1 - (clen - Math.min(slen, tlen))/(double)Math.max(slen,  tlen);
	}
	
	private static int[] combine(int[] s, int[] t) {
		int[] v = new int[s.length+t.length];
		for (int i = 0; i < s.length; i++) {
			v[i] = s[i];
		}
		for (int i = 0; i < t.length; i++) {
			v[s.length+i] = t[i];
		}
		return v;
	}
	
	private static int length(int[] values) {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		GZIPOutputStream zos;
		try {
			zos = new GZIPOutputStream(baos);

			for (int v: values) {
				zos.write(v);
			}
			zos.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return baos.size();
	}
}
