package org.rcsb.project3;

import java.util.List;

import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;

import scala.Tuple2;

/**
 * This class maps a pair of chains to the longest local common subsequence over the length of the chains
 * @author Chris Li
 */
public class SmithWatermanGotohP3 implements PairFunction<Tuple2<Integer,Integer>,String,Float> {

	private static final long serialVersionUID = 1L;
	private Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data = null;
    // print traceback if it is greater than 0
    private int traceback = 0;
    // open gap penalty
    private double open = 1;
    // extend gap penalty
    private double extend = 1;

    public SmithWatermanGotohP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data) {
		this.data = data;
	}
    
    public SmithWatermanGotohP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data, double open, double extend) {
		this.data = data;
		this.open = open;
		this.extend = extend;
	}
    
    /**
     * Constructor with traceback option
     * @param data
     * @param traceback
     */
	public SmithWatermanGotohP3(Broadcast<List<Tuple2<String,SequenceFeatureInterface<?>>>> data, int traceback) {
		this.data = data;
		this.traceback = traceback;
	}
	
	@Override
	public Tuple2<String, Float> call(Tuple2<Integer, Integer> tuple) {
		Tuple2<String,SequenceFeatureInterface<?>> t1 = this.data.getValue().get(tuple._1);
		Tuple2<String,SequenceFeatureInterface<?>> t2 = this.data.getValue().get(tuple._2);
		
		StringBuilder key = new StringBuilder();
		key.append(t1._1);
		key.append(",");
		key.append(t2._1);
		
		SequenceFeatureInterface<?> v1 = t1._2;
		SequenceFeatureInterface<?> v2 = t2._2;
		
		Alignment<?> SWAlignment = getAlignment(v1, v2, open, extend);
		float value = (float) SWAlignment.calculateScore();
		//value = value/v1.length();;
		//printTraceback(v1,v2,SWAlignment.getSequence1(),SWAlignment.getSequence2());
		int v1L = v1.length();
		int v2L = v2.length();
		if (v1L > v2L)
			value = value/v2L;
		else
			value = value/v1L;
		return new Tuple2<String, Float>(key.toString(), value);
	}
	
	@SuppressWarnings("unchecked")
	private <T,K> Alignment<T> getAlignment(SequenceFeatureInterface<T> v1,SequenceFeatureInterface<K> v2,double o, double e) {
		return SmithWatermanGotoh.align(v1, (SequenceFeatureInterface<T>)v2, o, e);
	}
	
	private void printTraceback(SequenceFeatureInterface<?> v1,SequenceFeatureInterface<?> v2,Integer[] v1Order,Integer[] v2Order) {
		//System.out.println("v1 length "+v1.length());
		//System.out.println("v2 length "+v2.length());
		//System.out.println("alignment length "+v1Order.length);
		System.out.print("v1 sequence ");
		for (int i = 0; i < v1Order.length; i++) {
			System.out.print(v1Order[i] + " ");
		}
		System.out.println();
		System.out.print("v2 sequence ");
		for (int i = 0; i < v2Order.length; i++) {
			System.out.print(v2Order[i] + " ");
		}
		System.out.println();
		System.out.println();
		/*int x = maxX;
		int y = maxY;
		String commonAngleV1 = " end at " + x;
		String commonAngleV2 = " end at " + y;
		while (x >= 0 && y >= 0 && b[x][y] >= 0) {
			if (b[x][y] == 0) {
				commonAngleV1 = v1.toString(x)+" \t"+commonAngleV1;
				commonAngleV2 = v2.toString(y)+" \t"+commonAngleV2;
				x--;
				y--;
			} else if (b[x][y] == 1) {
				commonAngleV1 = "xxx \t"+commonAngleV1;
				commonAngleV2 = "--- \t"+commonAngleV2;
				x--;
			} else if (b[x][y] == 2) {
				commonAngleV1 = "--- \t"+commonAngleV1;
				commonAngleV2 = "xxx \t"+commonAngleV2;
				y--;
			} else if (b[x][y] == -1) {
				break;
			}
		}
		if (x < 0)
			x++;
		if (y < 0)
			y++;
		commonAngleV1 = "start from "+ x + "\t" +commonAngleV1;
		commonAngleV2 = "start from "+ y + "\t" +commonAngleV2;
		System.out.println(commonAngleV1);
		System.out.println(commonAngleV2);*/
	}
}
