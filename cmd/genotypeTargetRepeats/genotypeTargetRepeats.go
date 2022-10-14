package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/realign"
	"github.com/guptarohit/asciigraph"
	"github.com/vertgenlab/gonomics/bed"
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
	"math"
	"os"
	"sort"
	"strings"
)

var debug int = 0

func usage() {
	fmt.Print(
		"genotypeTargetRepeats - Output a VCF of genotypes of targeted short simple repeats.\n\n" +
			"options:\n")
	flag.PrintDefaults()
}

// inputFiles is a custom type that gets filled by flag.Parse()
type inputFiles []string

// String to satisfy flag.Value interface
func (i *inputFiles) String() string {
	return strings.Join(*i, " ")
}

// Set to satisfy flag.Value interface
func (i *inputFiles) Set(value string) error {
	*i = append(*i, value)
	return nil
}

func main() {
	var inputs inputFiles
	flag.Var(&inputs, "i", "Input BAM file with alignments. Must be sorted and indexed. Can be declared more than once")
	var ref *string = flag.String("r", "", "Reference genome. Must be the same reference used for generating the BAM file.")
	var targets *string = flag.String("t", "", "BED file of targeted repeats. The 4th column must be the sequence of one repeat unit (e.g. CA for a CACACACA repeat), or 'RepeatLen'x'RepeatSeq' (e.g. 10xCA).")
	var output *string = flag.String("o", "stdout", "Output VCF file.")
	var bamOut *string = flag.String("bamOutPfx", "", "Output a BAM file with realigned reads. Only outputs reads that inform called genotypes. File will be name 'bamOutPfx'_'originalFilename'.")
	var targetPadding *int = flag.Int("tPad", 50, "Add INT bases of padding to either end of regions in targets file for selecting reads for realignment.")
	var minFlankOverlap *int = flag.Int("minFlank", 4, "A minimum of INT bases must be mapped on either side of the repeat to be considered an enclosing read.")
	var minMapQ *int = flag.Int("minMapQ", -1, "Minimum mapping quality (before realignment) to be considered for genotyping. Set to -1 for no filter.")
	var removeDups *bool = flag.Bool("removeDups", true, "Remove duplicate reads when genotyping.")
	var debugVal *int = flag.Int("debug", 0, "Set to 1 or greater for debug prints.")
	var minReads *int = flag.Int("minReads", 5, "Minimum total enclosing reads for genotyping.")
	flag.Parse()
	flag.Usage = usage

	if len(inputs) == 0 || *ref == "" {
		usage()
		log.Fatalln("ERROR: must input a VCF file with -i")
	}

	debug = *debugVal

	if *minMapQ > math.MaxUint8 {
		log.Fatalf("minMapQ out of range. max: %d\n", math.MaxUint8)
	}

	genotypeTargetRepeats(inputs, *ref, *targets, *output, *bamOut, *targetPadding, *minFlankOverlap, *minMapQ, *minReads, *removeDups)
}

func genotypeTargetRepeats(inputFiles []string, refFile, targetsFile, outputFile, bamOutPfx string, targetPadding, minFlankOverlap, minMapQ, minReads int, removeDups bool) {
	var err error
	ref := fasta.NewSeeker(refFile, "")
	targets := bed.Read(targetsFile)
	vcfOut := fileio.EasyCreate(outputFile)
	vcfHeader := generateVcfHeader(strings.Join(inputFiles, "\t"))
	vcf.NewWriteHeader(vcfOut, vcfHeader)

	// get bam reader for each file
	br := make([]*sam.BamReader, len(inputFiles))
	headers := make([]sam.Header, len(inputFiles))
	bamIdxs := make([]sam.Bai, len(inputFiles))
	for i := range inputFiles {
		br[i], headers[i] = sam.OpenBam(inputFiles[i])
		bamIdxs[i] = sam.ReadBai(inputFiles[i] + ".bai")
	}

	bamOutHandle := make([]io.WriteCloser, len(inputFiles))
	bamOut := make([]*sam.BamWriter, len(inputFiles))
	if bamOutPfx != "" {
		for i := range inputFiles {
			bamOutHandle[i] = fileio.EasyCreate(bamOutPfx + inputFiles[i])
			bamOut[i] = sam.NewBamWriter(bamOutHandle[i], headers[i])
		}
	}

	enclosingReads := make([][]*sam.Sam, len(inputFiles)) // first index is sample
	observedLengths := make([][]int, len(inputFiles))     // first index is sample
	var currVcf vcf.Vcf
	var i int
	alignerInput := make(chan sam.Sam)
	alignerOutput := realign.GoRealignIndels(alignerInput, ref)
	for _, region := range targets {
		for i = range inputFiles {
			enclosingReads[i], observedLengths[i] = getLenghtDist(enclosingReads[i], targetPadding, minMapQ, minFlankOverlap, removeDups, bamIdxs[i], region, br[i], bamOut[i], alignerInput, alignerOutput)
			if bamOutPfx != "" {
				for j := range enclosingReads {
					sam.WriteToBamFileHandle(bamOut[i], *enclosingReads[i][j], 0)
				}
			}
			slices.Sort(observedLengths[i])
		}

		if debug > 0 {
			fmt.Println(region, len(observedLengths), printLengths(observedLengths))
		}
		if debug > 1 {
			plot(observedLengths, minReads)
		}
		currVcf = callGenotypes(region, minReads, enclosingReads, observedLengths)
		vcf.WriteVcf(vcfOut, currVcf)
	}
	close(alignerInput)

	for i = range inputFiles {
		err = br[i].Close()
		exception.PanicOnErr(err)
		err = bamOut[i].Close()
		exception.PanicOnErr(err)
		err = bamOutHandle[i].Close()
	}
	err = ref.Close()
	exception.PanicOnErr(err)
	err = vcfOut.Close()
	exception.PanicOnErr(err)
}

func callGenotypes(region bed.Bed, minReads int, enclosingReads [][]*sam.Sam, observedLengths [][]int) vcf.Vcf {
	var ans vcf.Vcf

	return ans
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
		observedLengths[i] = calcRepeatLength(enclosingReads[i], region.ChromStart, region.ChromEnd, repeatSeq)
		if debug > 2 {
			fmt.Fprintln(os.Stderr, enclosingReads[i].QName, observedLengths[i], "start:", enclosingReads[i].Pos)
		}
	}
	return enclosingReads, observedLengths
}

func calcRepeatLength(read *sam.Sam, regionStart, regionEnd int, repeatSeq []dna.Base) int {
	var readIdx, refIdx, i int
	refIdx = int(read.Pos)

	// get to start of region
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
		repeatIdx--
		readIdx--
		refIdx--
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
	refIdx++
	// move forwards to calc repeat length
	var observedLength, maxLength int
	for refIdx < regionEnd && readIdx < len(read.Seq) {
		// move through repeat until mismatch
		for read.Seq[readIdx] == repeatSeq[repeatIdx] {
			observedLength++
			repeatIdx++
			readIdx++
			refIdx++
			if repeatIdx == len(repeatSeq) {
				repeatIdx = 0
			}
			if readIdx == len(read.Seq) {
				break
			}
		}
		if observedLength > maxLength {
			maxLength = observedLength
			observedLength = 0
		}
		// move forward until you get a base matching the repeat
		for readIdx < len(read.Seq) && read.Seq[readIdx] != repeatSeq[repeatIdx] {
			for repeatIdx = 0; repeatIdx < len(repeatSeq); repeatIdx++ {
				if read.Seq[readIdx] == repeatSeq[repeatIdx] {
					break
				}
			}
			if repeatIdx == len(repeatSeq) { // current read base does not match any base in repeat sequence
				repeatIdx = 0
				readIdx++
				refIdx++
			}
		}
	}
	return maxLength
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

func generateVcfHeader(samples string) vcf.Header {
	var header vcf.Header
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", samples))
	return header
}

func plot(observedLengths [][]int, minReads int) {
	readsPerSample := make([]int, len(observedLengths))
	p := make([][]float64, len(observedLengths))
	for i := range observedLengths {
		p[i] = make([]float64, 100)
		for j := range observedLengths[i] {
			p[i][observedLengths[i][j]]++
			readsPerSample[i]++
		}
	}
	for i := range p {
		//if readsPerSample[i] < minReads {
		//	continue
		//}
		fmt.Println(asciigraph.Plot(p[i], asciigraph.Height(5), asciigraph.Precision(0), asciigraph.SeriesColors(asciigraph.AnsiColor(i))))
	}
	fmt.Println(asciigraph.PlotMany(p, asciigraph.Precision(0), asciigraph.SeriesColors(
		asciigraph.Red,
		asciigraph.Yellow,
		asciigraph.Green,
		asciigraph.Blue,
		asciigraph.Cyan,
		asciigraph.BlueViolet,
		asciigraph.Brown,
		asciigraph.Gray,
		asciigraph.Orange,
		asciigraph.Olive,
	), asciigraph.Height(10)))
}

func printLengths(a [][]int) string {
	if len(a) == 0 {
		return ""
	}
	s := new(strings.Builder)
	for i := range a {
		if len(a[i]) == 0 {
			continue
		}
		s.WriteString(fmt.Sprintf("\t%d", a[i][0]))
		for j := 1; j < len(a[i]); j++ {
			s.WriteString(fmt.Sprintf(",%d", a[i][j]))
		}
	}
	return s.String()
}
