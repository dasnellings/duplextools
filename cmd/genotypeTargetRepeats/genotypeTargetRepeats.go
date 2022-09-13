package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/realign"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"math"
	"sort"
	"strings"
)

func usage() {
	fmt.Print(
		"genotypeTargetRepeats - Output a VCF of genotypes of targeted short simple repeats.\n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var input *string = flag.String("i", "", "Input BAM file with alignments. Must be sorted and indexed.")
	var ref *string = flag.String("r", "", "Reference genome. Must be the same reference used for generating the BAM file.")
	var targets *string = flag.String("t", "", "BED file of targeted repeats. The 4th column must be the sequence of one repeat unit (e.g. CA for a CACACACA repeat), or 'RepeatLen'x'RepeatSeq' (e.g. 10xCA).")
	var output *string = flag.String("o", "stdout", "Output VCF file.")
	var bamOut *string = flag.String("bamOut", "", "Output a BAM file with realigned reads. Only outputs reads that inform called genotypes.")
	var targetPadding *int = flag.Int("tPad", 50, "Add INT bases of padding to either end of regions in targets file for selecting reads for realignment.")
	var minFlankOverlap *int = flag.Int("minFlank", 4, "A minimum of INT bases must be mapped on either side of the repeat to be considered an enclosing read.")
	var minMapQ *int = flag.Int("minMapQ", -1, "Minimum mapping quality (before realignment) to be considered for genotyping. Set to -1 for no filter.")
	var removeDups *bool = flag.Bool("removeDups", true, "Remove duplicate reads when genotyping.")
	flag.Parse()
	flag.Usage = usage

	if *input == "" || *ref == "" {
		usage()
		log.Fatalln("ERROR: must input a VCF file with -i")
	}

	if *minMapQ > math.MaxUint8 {
		log.Fatalf("minMapQ out of range. max: %d\n", math.MaxUint8)
	}

	genotypeTargetRepeats(*input, *ref, *targets, *output, *bamOut, *targetPadding, *minFlankOverlap, *minMapQ, *removeDups)
}

func genotypeTargetRepeats(inputFile, refFile, targetsFile, outputFile, bamOutFile string, targetPadding, minFlankOverlap int, minMapQ int, removeDups bool) {
	var err error
	ref := fasta.NewSeeker(refFile, "")
	br, header := sam.OpenBam(inputFile)
	bamIdx := sam.ReadBai(inputFile + ".bai")
	targets := bed.Read(targetsFile)
	vcfOut := fileio.EasyCreate(outputFile)

	var bamOutHandle io.WriteCloser
	var bamOut *sam.BamWriter
	if bamOutFile != "" {
		bamOutHandle = fileio.EasyCreate(bamOutFile)
		bamOut = sam.NewBamWriter(bamOutHandle, header)
	}

	var enclosingReads []*sam.Sam
	var observedLengths []int
	alignerInput := make(chan sam.Sam)
	alignerOutput := realign.GoRealignIndels(alignerInput, ref)
	for _, region := range targets {
		enclosingReads, observedLengths = getLenghtDist(enclosingReads, targetPadding, minMapQ, minFlankOverlap, removeDups, bamIdx, region, br, bamOut, alignerInput, alignerOutput)
		fmt.Println(region, len(observedLengths), observedLengths)
	}
	close(alignerInput)

	err = ref.Close()
	exception.PanicOnErr(err)
	err = br.Close()
	exception.PanicOnErr(err)
	err = vcfOut.Close()
	exception.PanicOnErr(err)
	err = bamOut.Close()
	exception.PanicOnErr(err)
	err = bamOutHandle.Close()
}

func getLenghtDist(enclosingReads []*sam.Sam, targetPadding, minMapQ, minFlankOverlap int, removeDups bool, bamIdx sam.Bai, region bed.Bed, br *sam.BamReader, bamOut *sam.BamWriter, alignerInput chan<- sam.Sam, alignerOutput <-chan sam.Sam) ([]*sam.Sam, []int) {
	var start, end int
	var reads []sam.Sam
	enclosingReads = resetEnclosingReads(enclosingReads, len(reads)) // starts at len == 0, cap >= len(reads)

	// STEP 1: Find reads with initial alignment close to target as candidates for local realignment
	start = region.ChromStart - targetPadding
	end = region.ChromEnd + targetPadding
	if start < 0 {
		start = 0
	}
	reads = sam.SeekBamRegion(br, bamIdx, region.Chrom, uint32(start), uint32(end))
	if len(reads) == 0 {
		return enclosingReads, nil
	}

	// STEP 2: Realign reads to target region
	for i := range reads {
		if minMapQ != -1 && reads[i].MapQ < uint8(minMapQ) {
			continue
		}
		alignerInput <- reads[i]
		reads[i] = <-alignerOutput
	}

	// STEP 3: Determine which realigned reads overlap targets with the minimum flanking overlap
	for i := range reads {
		if minMapQ != -1 && reads[i].MapQ < uint8(minMapQ) {
			continue
		}
		if sam.IsUnmapped(reads[i]) {
			continue
		}
		if reads[i].GetChromStart() <= region.ChromStart-minFlankOverlap && reads[i].GetChromEnd() >= region.ChromEnd+minFlankOverlap {
			enclosingReads = append(enclosingReads, &reads[i])
			if bamOut != nil {
				sam.WriteToBamFileHandle(bamOut, reads[i], 0)
			}
		}
	}

	// STEP 4: Sort enclosing reads by position
	sort.Slice(enclosingReads, func(i, j int) bool {
		if enclosingReads[i].GetChromStart() < enclosingReads[j].GetChromStart() {
			return true
		}
		if enclosingReads[i].GetChromEnd() < enclosingReads[j].GetChromEnd() {
			return true
		}
		return true
	})

	// STEP 5: Remove duplicates
	if removeDups {
		enclosingReads = dedup(enclosingReads)
	}

	// STEP 6: Genotype repeats
	observedLengths := make([]int, len(enclosingReads))
	repeatSeq := parseRepeatSeq(region.Name)
	for i := range enclosingReads {
		observedLengths[i] = calcRepeatLength(enclosingReads[i], region.ChromStart, repeatSeq)
		fmt.Println(enclosingReads[i].QName, observedLengths[i], "start:", enclosingReads[i].Pos)
	}
	return enclosingReads, observedLengths
}

func calcRepeatLength(read *sam.Sam, regionStart int, repeatSeq []dna.Base) int {
	var readIdx, refIdx, i int
	refIdx = int(read.Pos)
	for i = range read.Cigar {
		if cigar.ConsumesReference(read.Cigar[i].Op) {
			refIdx += read.Cigar[i].RunLength
		}
		if cigar.ConsumesQuery(read.Cigar[i].Op) {
			readIdx += read.Cigar[i].RunLength
		}
		if refIdx >= regionStart {
			break
		}
	}
	if refIdx > regionStart {
		if cigar.ConsumesQuery(read.Cigar[i].Op) {
			readIdx -= refIdx - regionStart
		}
		refIdx -= refIdx - regionStart
	}
	readIdx++

	var repeatIdx int
	for repeatIdx = range repeatSeq {
		if read.Seq[readIdx] == repeatSeq[repeatIdx] {
			break
		}
	}

	// move backwards to look for misaligned repeat sequence
	for read.Seq[readIdx] == repeatSeq[repeatIdx] {
		repeatIdx -= 1
		readIdx -= 1
		if repeatIdx == -1 {
			repeatIdx = len(repeatSeq) - 1
		}
		if readIdx == -1 {
			break
		}
	}
	repeatIdx++
	if repeatIdx == len(repeatSeq) {
		repeatIdx = 0
	}
	readIdx++

	// move forwards to calc repeat length
	var observedLength int
	for read.Seq[readIdx] == repeatSeq[repeatIdx] {
		observedLength++
		repeatIdx += 1
		readIdx += 1
		if repeatIdx == len(repeatSeq) {
			repeatIdx = 0
		}
		if readIdx == len(read.Seq) {
			break
		}
	}

	return observedLength
}

func parseRepeatSeq(s string) []dna.Base {
	if strings.Contains(s, "x") {
		s = strings.Split(s, "x")[1]
	}
	return dna.StringToBases(s)
}

func dedup(reads []*sam.Sam) []*sam.Sam {
	for i := 1; i < len(reads); i++ {
		if reads[i].GetChromStart() == reads[i-1].GetChromStart() && reads[i].GetChromEnd() == reads[i-1].GetChromEnd() {
			slices.Delete(reads, i, i+1)
		}
	}
	return reads
}

func resetEnclosingReads(s []*sam.Sam, len int) []*sam.Sam {
	if cap(s) >= len {
		for i := range s {
			s[i] = nil
		}
		s = s[:0]
	} else {
		s = make([]*sam.Sam, 0, len)
	}
	return s
}
