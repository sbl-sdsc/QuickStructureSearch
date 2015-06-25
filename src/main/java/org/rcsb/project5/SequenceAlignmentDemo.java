package org.rcsb.project5;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SimpleSubstitutionMatrix;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;

public class SequenceAlignmentDemo {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ProteinSequence seq1  = null;
		try {
			seq1 = new ProteinSequence("MSTNPKPQRKTKRNTNRRPQDVKFPGG");
		} catch (CompoundNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		ProteinSequence seq2  = null;
		try {
//			seq2 = new ProteinSequence("PKPQRKTKRNTNRRPQDVK");
			seq2 = new ProteinSequence("PKPQRKTKRNPNRRPQDVK");
		} catch (CompoundNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
        SubstitutionMatrix<AminoAcidCompound> matrix = SimpleSubstitutionMatrix.getBlosum62();

        SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(seq1, seq2,
                PairwiseSequenceAlignerType.GLOBAL, new SimpleGapPenalty(), matrix);
        System.out.printf("%n%s vs %s%n%s", pair.getQuery().getAccession(), pair.getTarget().getAccession(), pair);
        System.out.println("query: " + pair.getCompoundInQueryAt(1));
        System.out.println("target: " + pair.getCompoundInTargetAt(1));
	}
}
