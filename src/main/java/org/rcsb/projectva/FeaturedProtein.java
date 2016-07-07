package org.rcsb.projectva;
import java.util.Map;


import java.util.Map.Entry;
import java.io.ByteArrayOutputStream;
import java.util.zip.GZIPOutputStream;

import org.rcsb.project3.*;


public class FeaturedProtein {
	/*
	public static void main(String[] args) {
		if(args.length > 0) {
			m_seqInfo = args[0];
		}
		
	}
	*/
	public static<T> int numberFeaturedProtein(T feature)
	{
		System.out.println("sequence info = " + m_seqInfo);
		
		return 0;
	}
	
	public static<T> boolean ifProteinFeatured(T feature, SequenceFeatureInterface<T> seq)
	{
		boolean result = false;
		for (int i = 0; i < seq.length(); i++) {
			if (seq.get(i).equals(feature)) result = true;
		}
		
		return result;
			
	}
	
	static String m_seqInfo;
}
