package realign

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

var gapOpen int64 = -600
var gapExtend int64 = -20

func GoRealignIndels(reads <-chan sam.Sam, ref *fasta.Seeker) <-chan sam.Sam {
	output := make(chan sam.Sam, 1000)
	go realignIndelsEngine(reads, output, ref)
	return output
}

func realignIndels(in <-chan sam.Sam, out chan<- sam.Sam, ref *fasta.Seeker) {
	var currStart, currEnd int
	var currRegion []dna.Base
	var score int64
	var cig []align.Cigar

	for r := range in {
		if !(r.GetChromStart() >= currStart+200 && r.GetChromEnd() <= currEnd-200) {
			currStart, currEnd, currRegion = getRegion(r, ref)
			dna.AllToUpper(currRegion)
		}
		score, cig = align.AffineGapLocal(currRegion, r.Seq, align.HumanChimpTwoScoreMatrix, gapOpen, gapExtend)
		updateRead(&r, cig, currStart, currEnd, score)
		out <- r
	}
	close(out)
}

func realignIndelsEngine(in <-chan sam.Sam, out chan<- sam.Sam, ref *fasta.Seeker) {
	var currStart, currEnd int
	var currRegion []dna.Base
	var packet align.TargetQueryPair
	inputs, outputs := align.GoAffineGapLocalEngine(align.HumanChimpTwoScoreMatrix, gapOpen, gapExtend)

	for r := range in {
		if !(r.GetChromStart() >= currStart+200 && r.GetChromEnd() <= currEnd-200) {
			currStart, currEnd, currRegion = getRegion(r, ref)
			dna.AllToUpper(currRegion)
		}

		packet.Target = currRegion
		packet.Query = r.Seq
		inputs <- packet
		packet = <-outputs
		updateRead(&r, packet.Cigar, currStart, currEnd, packet.Score)
		out <- r
	}
	close(out)
}

func getRegion(read sam.Sam, ref *fasta.Seeker) (start, end int, region []dna.Base) {
	var err error
	var pad int = 1000
	start = read.GetChromStart() - pad
	end = read.GetChromEnd() + pad
	if start < 0 {
		start = 0
	}
	region, err = fasta.SeekByName(ref, read.RName, start, end)
	if err != nil {
		log.Printf("ERROR: problem fetching reference %s:%d-%d for read %s\n", read.RName, start, end, read.QName)
		log.Fatalln(err)
	}

	return
}

func cigConv(c []align.Cigar) []cigar.Cigar {
	ans := make([]cigar.Cigar, len(c))
	for i := range c {
		switch c[i].Op {
		case align.ColM:
			ans[i].Op = 'M'
		case align.ColI:
			ans[i].Op = 'I'
		case align.ColD:
			ans[i].Op = 'D'
		}
		ans[i].RunLength = int(c[i].RunLength)
	}
	return ans
}

func updateRead(r *sam.Sam, cig []align.Cigar, cigStart, cigEnd int, score int64) {
	var alignStart int
	alignStart = cigStart
	if cig[0].Op == align.ColD {
		alignStart += int(cig[0].RunLength)
		cig = cig[1:]
	}
	if cig[len(cig)-1].Op == align.ColD {
		cig = cig[:len(cig)-1]
	}
	r.Pos = uint32(alignStart) + 1
	r.Cigar = cigConv(cig)
	//r.Extra += fmt.Sprintf("\tSC:i:%d", score)
}
