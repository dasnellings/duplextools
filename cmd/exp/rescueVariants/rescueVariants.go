package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/barcode"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"path"
	"sort"
	"strings"
)

func usage() {
	fmt.Print(
		"rescueVariants - Re-evaluate the presence of high-confidence variants in read families with less support.\n\n" +
			"Usage:\n" +
			"  rescueVariants [options] -v highConfidence.vcf -b alignedReads.bam -r reference.fasta > output.vcf\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func main() {
	ref := flag.String("r", "", "`FASTA` file with reference genome used to align input BAM. Must be indexed.")
	highConfFile := flag.String("v", "", "`VCF` file containing high confidence dsSNVs to attempt rescue.")
	bamFile := flag.String("b", "", "`BAM` file to attempt rescue from.")
	output := flag.String("o", "stdout", "Output `VCF` file containing all high confidence dsSNVs, in addition to any rescued variants.")
	minReads := flag.Int("minReads", 2, "Minimum reads for read family to be considered positive for variant.")
	endPad := flag.Int("ignoreEnds", 3, "Ignore bases within # of end of a read.")
	minMapQ := flag.Int("minMapQ", 20, "Minimum mapping quality.")
	maxSoftClipFraction := flag.Float64("maxSoftClipFraction", 0.2, "Maximum fraction of read that may be soft clipped.")
	countOverlappingPairs := flag.Bool("countOverlappingPairs", false, "Count both reads in overlapping regions of read pairs. By only 1 base is contributed in overlapping regions of read pairs.")
	allowSuppAln := flag.Bool("allowSupplementaryAlignments", false, "Allow variants using reads that have supplementary alignments annotated.")
	includeSelf := flag.Bool("includeSelf", false, "When ignoreSelf is false, reads from the read family the variant was discovered in will be ignored.")
	minAf := flag.Float64("minAF", 0.9, "Minimum fraction of reads with alternate allele **Within a read family and within strand** to be considered a variant.")
	minBaseQuality := flag.Int("minBaseQuality", 30, "Minimum base quality to be considered for calling. Bases below threshold will be ignored.")
	baseQualPenalty := flag.Float64("baseQualPenalty", 0.5, "Penalty for positions with low quality base. Each read with a base < minBaseQuality counts towards baseQualPenalty fraction of a read for allele frequency calculations. Note that low quality bases are N-masked and so will always count AGAINST the alternate allele. (e.g. by default each read with a low quality base counts as 0.5 reads for allele frequency determination.")
	excludeRescued := flag.Bool("excludeRescued", false, "Do not output any sites with >1 rescued variants.")
	excludeUnrescued := flag.Bool("excludeUnrescued", false, "Do not output any sites with 0 rescued variants.")
	flag.Parse()
	flag.Usage = usage

	if *ref == "" || *highConfFile == "" || *bamFile == "" {
		flag.Usage()
		log.Fatalln("ERROR: Must input values for -v, -b, and -r")
	}

	rescueVariants(*ref, *highConfFile, *bamFile, *output, *endPad, *minMapQ, *minReads, *countOverlappingPairs, *allowSuppAln, *includeSelf, *minAf, *minBaseQuality, *baseQualPenalty, *maxSoftClipFraction, *excludeRescued, *excludeUnrescued)
}

func rescueVariants(refFile, highConfFile, bamFile, outFile string, endPad, minMapQ, minReads int, countOverlappingPairs, allowSuppAln, includeSelf bool, minAf float64, minBaseQuality int, baseQualPenalty, maxSoftClipFraction float64, excludeRescued, excludeUnrescued bool) {
	ref := fasta.NewSeeker(refFile, "")
	defer cleanup(ref)
	variants, vcfHeader := vcf.GoReadToChan(highConfFile)
	if !excludeRescued {
		vcfHeader = updateHeader(vcfHeader)
	}
	out := fileio.EasyCreate(outFile)
	defer cleanup(out)
	vcf.NewWriteHeader(out, vcfHeader)
	bam, bamHeader := sam.OpenBam(bamFile)
	bai := sam.ReadBai(bamFile + ".bai")
	defer cleanup(bam)

	var reads []sam.Sam
	var familyReads [][]sam.Sam
	var familyIds []string
	var pile sam.Pile
	var variantType, originFamilyId string
	var refSeq, altSeq []dna.Base
	var i, depth, altCount, positiveFamilyCount, totalRescuedSites, sitesRemoved int
	for v := range variants {
		reads = sam.SeekBamRegionRecycle(bam, bai, v.Chr, uint32(v.Pos-1), uint32(v.Pos), reads)
		originFamilyId = getFamily(v)

		refSeq = dna.StringToBases(v.Ref)
		altSeq = dna.StringToBases(v.Alt[0])

		switch {
		case len(refSeq) == len(altSeq):
			variantType = "SNV"
		case len(refSeq) > len(altSeq):
			variantType = "Deletion"
		case len(refSeq) < len(altSeq):
			variantType = "Insertion"
		}

		positiveFamilyCount = 0
		familyReads, familyIds = splitReadsByFamily(reads, familyReads, familyIds, allowSuppAln, endPad, minBaseQuality, minMapQ, originFamilyId, includeSelf, maxSoftClipFraction)
		for i = range familyReads {
			pile = getPile(familyReads[i], bamHeader, v.Pos, variantType, countOverlappingPairs)
			depth = calcDepth(pile, baseQualPenalty)

			switch variantType {
			case "SNV":
				altCount = pile.CountF[altSeq[0]] + pile.CountR[altSeq[0]]
			case "Deletion":
				altCount = pile.DelCountF[len(refSeq)-len(altSeq)] + pile.DelCountR[len(refSeq)-len(altSeq)]
			case "Insertion":
				altCount = pile.InsCountF[v.Alt[0][1:]] + pile.InsCountR[v.Alt[0][1:]] // need to remove leading base for key to map
			}
			if altCount < minReads || float64(altCount)/float64(depth) < minAf {
				continue
			}

			positiveFamilyCount++
		}
		if !excludeRescued && v.Info != "" {
			v.Info += ";"
		}

		if positiveFamilyCount == 0 && excludeUnrescued {
			sitesRemoved++
			continue
		}
		if positiveFamilyCount > 0 && excludeRescued {
			sitesRemoved++
			continue
		}
		if !excludeRescued {
			v.Info += fmt.Sprintf("RescueFamilyCount=%d", positiveFamilyCount)
		}
		vcf.WriteVcf(out, v)

		if positiveFamilyCount > 0 {
			totalRescuedSites++
		}
	}
	if !excludeRescued {
		log.Printf("Total Mosaic Sites: %d\t%s\t%s\n", totalRescuedSites, path.Base(highConfFile), path.Base(bamFile))
	} else if excludeRescued || excludeUnrescued {
		log.Printf("Sites Removed: %d\n", sitesRemoved)
	}
}

func splitReadsByFamily(reads []sam.Sam, familyReads [][]sam.Sam, familyIds []string, allowSuppAln bool, endPad, minBaseQuality, minMapQ int, originFamilyId string, includeSelf bool, maxSoftClipFraction float64) ([][]sam.Sam, []string) {
	familyIds = familyIds[:0]
	for i := range familyReads {
		familyReads[i] = familyReads[i][:0]
	}
	familyReads = familyReads[:0]

	var famId string
	var j int
	var foundFam bool
	for i := range reads {
		if reads[i].MapQ < uint8(minMapQ) || (hasSuppAln(reads[i]) && !allowSuppAln) {
			continue
		}
		if softClipFraction(&reads[i]) > maxSoftClipFraction {
			continue
		}
		clipReadEnds(&reads[i], endPad)
		maskLowQualityBases(&reads[i], minBaseQuality)

		sam.ParseExtra(&reads[i])
		famId = barcode.GetRF(&reads[i])

		if !includeSelf && famId == originFamilyId {
			continue
		}

		foundFam = false
		for j = range familyIds {
			if famId == familyIds[j] {
				foundFam = true
				break
			}
		}
		if foundFam {
			familyReads[j] = append(familyReads[j], reads[i])
		} else {
			familyIds = append(familyIds, famId)
			familyReads = append(familyReads, []sam.Sam{reads[i]})
		}
	}

	for k := range familyReads {
		sort.Slice(familyReads[k], func(i, j int) bool {
			return familyReads[k][i].Pos < familyReads[k][j].Pos
		})
	}
	return familyReads, familyIds
}

func getPile(reads []sam.Sam, header sam.Header, pos int, variantType string, countOverlappingPairs bool) sam.Pile {
	if len(reads) == 0 {
		return sam.Pile{}
	}

	samChan := make(chan sam.Sam, len(reads))
	for i := range reads {
		sclipTerminalIns(&reads[i])
		samChan <- reads[i]
	}
	close(samChan)

	pileChan := sam.GoPileup(samChan, header, false, nil, nil)
	for p := range pileChan {
		if (variantType == "SNV" || variantType == "Insertion") && p.Pos == uint32(pos) {
			if !countOverlappingPairs {
				removeBasesFromOverlappingReadPairs(&p)
			}
			return p
		}
		if variantType == "Deletion" && p.Pos == uint32(pos+1) {
			return p
		}
	}
	return sam.Pile{}
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}

// calcDepth returns the number of reads in the input pile
func calcDepth(p sam.Pile, baseQualPenalty float64) int {
	var depth int
	var maskCount int
	for i := range p.CountF {
		if i == int(dna.N) {
			maskCount += p.CountF[i] + p.CountR[i]
			continue
		}
		depth += p.CountF[i] + p.CountR[i]
	}
	depth += int(float64(maskCount) * baseQualPenalty)
	return depth
}

func updateHeader(h vcf.Header) vcf.Header {
	newLine := "##INFO=<ID=RescueFamilyCount,Number=1,Type=Integer,Description=\"Number of read families harboring mutation using relaxed calling from the duplexTools/rescueVariants command.\">"
	for i := range h.Text {
		if strings.HasPrefix(h.Text[i], "##INFO") {
			h.Text = slices.Insert(h.Text, i, newLine)
			return h
		}
	}
	// no info lines already present
	h.Text = slices.Insert(h.Text, len(h.Text)-1, newLine)
	return h
}

func clipReadEnds(s *sam.Sam, clipLen int) {
	if s.Cigar == nil || len(s.Cigar) == 0 || s.Cigar[0].Op == '*' {
		return
	}

	var anyNonClip bool
	for i := range s.Cigar {
		if s.Cigar[i].Op != 'S' {
			anyNonClip = true
			break
		}
	}

	if !anyNonClip {
		return
	}

	clipFwd(s, clipLen)
	clipRev(s, clipLen)

	// collapse cigar if everything is soft clipped
	if len(s.Cigar) == 2 && s.Cigar[0].Op == 'S' && s.Cigar[1].Op == 'S' {
		s.Cigar[0].RunLength += s.Cigar[1].RunLength
		s.Cigar = s.Cigar[:1]
	}

	//if cigar.QueryLength(s.Cigar) != len(s.Seq) {
	//	log.Panic("something went horribly wrong with cigar\n", s)
	//}
}

func clipFwd(s *sam.Sam, clipLen int) {
	if clipLen < 1 {
		return
	}

	// check if first index is soft clip, if not make a soft clip with len = 0
	if s.Cigar[0].Op != 'S' {
		s.Cigar = slices.Insert(s.Cigar, 0, cigar.Cigar{Op: 'S', RunLength: 0})
	}
	var numToClip int = clipLen
	var currNumToClip int
	for i := 1; numToClip > 0; i++ {
		// increment pos as well as cigar
		switch s.Cigar[i].Op {
		case 'M':
			currNumToClip = min(s.Cigar[i].RunLength, numToClip)
			s.Cigar[i].RunLength -= currNumToClip
			s.Cigar[0].RunLength += currNumToClip
			s.Pos += uint32(currNumToClip)
			numToClip -= currNumToClip

		case 'D':
			s.Pos += uint32(s.Cigar[i].RunLength)
			s.Cigar[i].RunLength = 0

		case 'I':
			currNumToClip = min(s.Cigar[i].RunLength, numToClip)
			s.Cigar[0].RunLength += currNumToClip
			s.Cigar[i].RunLength -= currNumToClip
			numToClip -= currNumToClip

		case 'S':
			s.Cigar = cleanCigar(s.Cigar)
			return
		}
	}
	s.Cigar = cleanCigar(s.Cigar)
}

func clipRev(s *sam.Sam, clipLen int) {
	if clipLen < 1 {
		return
	}

	// check if last index is soft clip, if not make a soft clip with len = 0
	if s.Cigar[len(s.Cigar)-1].Op != 'S' {
		s.Cigar = append(s.Cigar, cigar.Cigar{Op: 'S', RunLength: 0})
	}
	var numToClip int = clipLen
	var currNumToClip int
	lastIdx := len(s.Cigar) - 1
	for i := lastIdx - 1; numToClip > 0; i-- {
		// increment pos as well as cigar
		switch s.Cigar[i].Op {
		case 'M', 'I':
			currNumToClip = min(s.Cigar[i].RunLength, numToClip)
			s.Cigar[i].RunLength -= currNumToClip
			s.Cigar[lastIdx].RunLength += currNumToClip
			numToClip -= currNumToClip

		case 'D':
			s.Cigar[i].RunLength = 0

		case 'S':
			s.Cigar = cleanCigar(s.Cigar)
			return
		}
	}
	s.Cigar = cleanCigar(s.Cigar)
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

func cleanCigar(c []cigar.Cigar) []cigar.Cigar {
	// remove all indexes with RunLength of 0
	for i := 0; i < len(c); i++ {
		if c[i].RunLength == 0 {
			c = slices.Delete(c, i, i+1)
			i--
		}
	}
	return c
}

func hasSuppAln(r sam.Sam) bool {
	_, found, err := sam.QueryTag(r, "SA")
	if err != nil || !found {
		return false
	}
	return true
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func removeBasesFromOverlappingReadPairs(p *sam.Pile) {
	for i := range p.CountF {
		if p.CountF[i] > p.CountR[i] {
			p.CountR[i] = 0
		} else {
			p.CountF[i] = 0
		}
	}

	for key := range p.DelCountF {
		if p.DelCountF[key] > p.DelCountR[key] {
			p.DelCountR[key] = 0
		} else {
			p.DelCountF[key] = 0
		}
	}

	for key := range p.DelCountR {
		if p.DelCountF[key] > p.DelCountR[key] {
			p.DelCountR[key] = 0
		} else {
			p.DelCountF[key] = 0
		}
	}

	for key := range p.InsCountF {
		if p.InsCountF[key] > p.InsCountR[key] {
			p.InsCountR[key] = 0
		} else {
			p.InsCountF[key] = 0
		}
	}

	for key := range p.InsCountR {
		if p.InsCountF[key] > p.InsCountR[key] {
			p.InsCountR[key] = 0
		} else {
			p.InsCountF[key] = 0
		}
	}
}

func sclipTerminalIns(s *sam.Sam) {
	if len(s.Cigar) == 0 || s.Cigar[0].Op == '*' {
		return
	}
	if s.Cigar[0].Op == 'I' {
		s.Cigar[0].Op = 'S'
	}
	if s.Cigar[len(s.Cigar)-1].Op == 'I' {
		s.Cigar[len(s.Cigar)-1].Op = 'S'
	}

	// catch case where beginning/end of read is already soft clipped
	if len(s.Cigar) >= 2 && s.Cigar[0].Op == 'S' && s.Cigar[1].Op == 'I' {
		s.Cigar[1].Op = 'S'
		s.Cigar[1].RunLength += s.Cigar[0].RunLength
		s.Cigar = s.Cigar[1:]
	}

	if len(s.Cigar) >= 2 && s.Cigar[len(s.Cigar)-1].Op == 'S' && s.Cigar[len(s.Cigar)-2].Op == 'I' {
		s.Cigar[len(s.Cigar)-2].Op = 'S'
		s.Cigar[len(s.Cigar)-2].RunLength += s.Cigar[len(s.Cigar)-1].RunLength
		s.Cigar = s.Cigar[:len(s.Cigar)-1]
	}
}

func getFamily(v vcf.Vcf) string {
	idx := -1
	for i := range v.Format {
		if v.Format[i] == "RF" {
			idx = i
			break
		}
	}

	if idx == -1 {
		return ""
	}

	return v.Samples[0].FormatData[idx]
}

func softClipFraction(r *sam.Sam) float64 {
	totalLen := len(r.Seq)
	var sClipCount int
	for i := range r.Cigar {
		if r.Cigar[i].Op == 'S' {
			sClipCount += r.Cigar[i].RunLength
		}
	}
	return float64(sClipCount) / float64(totalLen)
}
