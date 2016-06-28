package org.rcsb.project10;

import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFunction;

import scala.Tuple2;

/**
 * This class takes a tuple of scores and returns a binary classification.
 * The constructor takes a threshold value between [0,1] that use use to
 * decide if a score is considered negative (< threshold) or positive (>= threshold).
 * TP: true positive, FP: false positive, TN: true negative, FN: false negative
 * For details see: https://en.wikipedia.org/wiki/Precision_and_recall
 * 
 * @author Peter Rose
 *
 */
public class BinaryClassificationMapper  implements PairFunction<Tuple2<String, Tuple2<Double,Double>>, String, String> {
	private static final long serialVersionUID = 1847235980329917612L;
	private double threshold;
    
	public BinaryClassificationMapper(double threshold) {
		this.threshold = threshold;
	}
	
	@Override
	public Tuple2<String, String> call(Tuple2<String, Tuple2<Double, Double>> t) {
	    double v1 = t._2._1; // first value in tuple
	    double v2 = t._2._2; // second value in tuple
	
		String classification = "";
        if (v1 >= this.threshold) {
        	if (v2 >= this.threshold) {
        		classification = "TP";
        	} else {
        		classification = "FN";
        	}
        } else if (v1 < this.threshold) {
        	if (v2 < this.threshold) {
        		classification = "TN";
        	} else {
                classification = "FP";
        	}
        }
        
        return new Tuple2<String,String>(t._1, classification);
	}
}
