package repeats

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
)

func FindPerfectRepeat(reference *fasta.Seeker, r *Record) (bed.Bed, int, []dna.Base) {
	refseq, err := fasta.SeekByName(reference, r.Chr, r.Start, r.End)
	exception.PanicOnErr(err)
	dna.AllToUpper(refseq)

	bestStart, bestEnd, numRepeats := bestMatch(refseq, r.Seq)
	var ans bed.Bed
	ans.FieldsInitialized = 4
	ans.Chrom = r.Chr
	ans.ChromStart = r.Start + bestStart
	ans.ChromEnd = r.Start + bestEnd
	ans.Name = fmt.Sprintf("%dx%s", numRepeats, dna.BasesToString(r.Seq))
	return ans, numRepeats, r.Seq
}

func incrementPatternIdx(pattern []dna.Base, idx *int) {
	if *idx == len(pattern)-1 {
		*idx = 0
		return
	}
	*idx++
}

func bestMatch(seq, pattern []dna.Base) (bestStart, bestEnd, numRepeats int) {
	var currMatchStart, bestMatchStart, bestMatchEnd, currMatchEnd, currPatternIdx, lenMatching, i int
	for i = 0; i < len(seq); i++ {
		if seq[i] == pattern[currPatternIdx] {
			if lenMatching == 0 {
				currMatchStart = i
			}
			lenMatching++
			incrementPatternIdx(pattern, &currPatternIdx)
		} else {
			currMatchEnd = i - currPatternIdx
			if bestMatchEnd-bestMatchStart < currMatchEnd-currMatchStart {
				bestMatchStart = currMatchStart
				bestMatchEnd = currMatchEnd
			}
			currPatternIdx = 0
			lenMatching = 0
			currMatchStart = i
			if seq[i] == pattern[currPatternIdx] {
				incrementPatternIdx(pattern, &currPatternIdx)
			}
		}
	}

	if bestMatchEnd-bestMatchStart < i-currMatchStart {
		bestMatchStart = currMatchStart
		bestMatchEnd = i
	}

	return bestMatchStart, bestMatchEnd, (bestMatchEnd - bestMatchStart) / len(pattern)
}
