package filter

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"sort"
)

func Filter(input, output, genomicBam, genomicVcf, snpVcf string, minCoverage int, maxReadFrac float64, maxReads int, minBaseQuality int) {
	var err error
	out := fileio.EasyCreate(output)
	inChan, header := vcf.GoReadToChan(input)
	vcf.NewWriteHeader(out, header)
	gBam, gBamHeader := sam.OpenBam(genomicBam)
	gBai := sam.ReadBai(genomicBam + ".bai")

	var gVcfChan <-chan vcf.Vcf
	if genomicVcf != "" {
		gVcfChan, _ = vcf.GoReadToChan(genomicVcf)
	}

	var excludeIntervals []interval.Interval
	var padding int
	if genomicVcf != "" {
		for v := range gVcfChan {
			padding = numbers.Max(len(v.Ref), len(v.Alt[0]))
			excludeIntervals = append(excludeIntervals, bed.Bed{Chrom: v.Chr, ChromStart: v.Pos - padding, ChromEnd: v.Pos + padding, FieldsInitialized: 3})
		}
	}

	if snpVcf != "" {
		snpVcfChan, _ := vcf.GoReadToChan(snpVcf)
		for v := range snpVcfChan {
			padding = numbers.Max(len(v.Ref), len(v.Alt[0]))
			excludeIntervals = append(excludeIntervals, bed.Bed{Chrom: v.Chr, ChromStart: v.Pos - padding, ChromEnd: v.Pos + padding, FieldsInitialized: 3})
		}
	}

	var tree map[string]*interval.IntervalNode
	if genomicVcf != "" {
		tree = interval.BuildTree(excludeIntervals)
	}

	filterGermline(inChan, out, gBam, gBamHeader, gBai, tree, minCoverage, maxReadFrac, maxReads, minBaseQuality)

	err = gBam.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func filterGermline(inChan <-chan vcf.Vcf, out *fileio.EasyWriter, gBam *sam.BamReader, gBamHeader sam.Header, gBai sam.Bai, excludeTree map[string]*interval.IntervalNode, minCoverage int, maxReadFrac float64, maxReadsLimit int, minBaseQuality int) {
	var reads []sam.Sam
	var p sam.Pile
	var maxReads, obsReads, delLen int
	var altBase dna.Base
	var prevChr string
	// loop over each variant in input and determine if it belongs in output
	for v := range inChan {
		if prevChr == "" {
			prevChr = v.Chr
		}

		if v.Chr != prevChr {
			if !checkChr(v.Chr, gBamHeader.Chroms) {
				log.Printf("WARNING: %s present in VCF, but not in bam file. Skipping variant.\n", v.Chr)
				continue
			}
			prevChr = v.Chr
		}

		if excludeTree != nil && len(interval.Query(excludeTree, v, "any")) > 0 {
			continue
		}

		p, reads = retrievePile(v, gBam, gBai, gBamHeader, reads, minBaseQuality)
		log.Printf("running %s\t%d\t%s\t%s\tReads:%d\tA:%d\tC:%d\tG:%d\tT:%d\tGap:%d\tIns:%d\tDel:%d\n", v.Chr, v.Pos, v.Ref, v.Alt[0], len(reads),
			p.CountF[dna.A]+p.CountR[dna.A],
			p.CountF[dna.C]+p.CountR[dna.C],
			p.CountF[dna.G]+p.CountR[dna.G],
			p.CountF[dna.T]+p.CountR[dna.T],
			p.CountF[dna.Gap]+p.CountR[dna.Gap],
			sumIns(p),
			sumDel(p))
		if len(reads) < minCoverage {
			continue
		}
		maxReads = int(float64(len(reads)) * maxReadFrac)
		if maxReads < 1 {
			maxReads = 1
		}

		// determine variant type and get number of reads in reference
		switch {
		case len(v.Ref) == 1 && len(v.Alt[0]) == 1: // substitution
			altBase = dna.StringToBase(v.Alt[0])
			obsReads = p.CountF[altBase] + p.CountR[altBase]

		case len(v.Ref) > len(v.Alt[0]): // deletion
			delLen = len(v.Ref) - len(v.Alt[0])
			obsReads = p.DelCountF[delLen] + p.DelCountR[delLen]

		case len(v.Alt[0]) > len(v.Ref): // insertion
			obsReads = p.InsCountF[v.Alt[0][1:]] + p.InsCountR[v.Alt[0][1:]]

		default:
			log.Panicf("something went horribly wrong with the following variant\t%s\n", v)
		}

		if obsReads > maxReads || obsReads > maxReadsLimit {
			continue
		}

		vcf.WriteVcf(out, v)
	}
}

func retrievePile(v vcf.Vcf, gBam *sam.BamReader, gBai sam.Bai, gBamHeader sam.Header, reads []sam.Sam, minBaseQuality int) (sam.Pile, []sam.Sam) {
	start := uint32(v.Pos) - 1
	stop := uint32(v.Pos)
	pos := v.Pos
	if len(v.Ref) > 1 { // deletion
		start++
		stop++
		pos++
	}

	reads = sam.SeekBamRegionRecycle(gBam, gBai, v.Chr, start, stop, reads)
	for i := range reads {
		maskLowQualityBases(&reads[i], minBaseQuality)
	}
	sort.Slice(reads, func(i, j int) bool { return reads[i].Pos < reads[j].Pos })
	piles := filterPileup(reads, gBamHeader)
	for i := range piles {
		if int(piles[i].Pos) == pos {
			return piles[i], reads
		}
	}
	return sam.Pile{}, reads
}

func filterPileup(reads []sam.Sam, header sam.Header) []sam.Pile {
	if len(reads) == 0 {
		return nil
	}

	samChan := make(chan sam.Sam, len(reads))
	for i := range reads {
		samChan <- reads[i]
	}
	close(samChan)

	ans := make([]sam.Pile, 0, 100)
	// TODO terribly inefficient to get piles for the whole region when we could smartly get the individual pile, but it's fast enough for now
	pileChan := sam.GoPileup(samChan, header, false, nil, nil)
	for p := range pileChan {
		ans = append(ans, p)
	}
	return ans
}

func sumIns(p sam.Pile) int {
	var ans int
	for _, val := range p.InsCountF {
		ans += val
	}
	for _, val := range p.InsCountR {
		ans += val
	}
	return ans
}

func sumDel(p sam.Pile) int {
	var ans int
	for _, val := range p.DelCountF {
		ans += val
	}
	for _, val := range p.DelCountR {
		ans += val
	}
	return ans
}

func checkChr(chr string, list []chromInfo.ChromInfo) bool {
	for i := range list {
		if list[i].Name == chr {
			return true
		}
	}
	return false
}

func maskLowQualityBases(s *sam.Sam, minQual int) {
	var currQual uint8
	for i := range s.Qual {
		currQual = s.Qual[i] - 33
		if currQual < uint8(minQual) {
			s.Seq[i] = dna.N
		}
	}
}
