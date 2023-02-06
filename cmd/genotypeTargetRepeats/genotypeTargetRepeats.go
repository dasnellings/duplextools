package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/gmm"
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
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
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
	var allowDups *bool = flag.Bool("allowDups", false, "Do not remove duplicate reads when genotyping.")
	var debugVal *int = flag.Int("debug", 0, "Set to 1 or greater for debug prints.")
	var minReads *int = flag.Int("minReads", 5, "Minimum total enclosing reads for genotyping.")
	var alignerThreads *int = flag.Int("alnThreads", 1, "Number of alignment threads.")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to `file`")
	memprofile := flag.String("memprofile", "", "write memory profile to `file`")
	flag.Parse()
	flag.Usage = usage

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	if len(inputs) == 0 || *ref == "" {
		usage()
		log.Fatalln("ERROR: must input a VCF file with -i")
	}

	debug = *debugVal

	if *minMapQ > math.MaxUint8 {
		log.Fatalf("minMapQ out of range. max: %d\n", math.MaxUint8)
	}

	genotypeTargetRepeats(inputs, *ref, *targets, *output, *bamOut, *targetPadding, *minFlankOverlap, *minMapQ, *minReads, !*allowDups, *alignerThreads)

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}

func genotypeTargetRepeats(inputFiles []string, refFile, targetsFile, outputFile, bamOutPfx string, targetPadding, minFlankOverlap, minMapQ, minReads int, removeDups bool, alignerThreads int) {
	var err error
	var ref *fasta.Seeker
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
		if _, err = os.Stat(inputFiles[i] + ".bai"); !errors.Is(err, os.ErrNotExist) {
			bamIdxs[i] = sam.ReadBai(inputFiles[i] + ".bai")
		} else {
			bamIdxs[i] = sam.ReadBai(strings.TrimSuffix(inputFiles[i], ".bam") + ".bai")
		}
	}

	bamOutHandle := make([]io.WriteCloser, len(inputFiles))
	bamOut := make([]*sam.BamWriter, len(inputFiles))
	if bamOutPfx != "" {
		for i := range inputFiles {
			words := strings.Split(inputFiles[i], "/")
			words[len(words)-1] = bamOutPfx + "_" + words[len(words)-1]
			bamOutHandle[i] = fileio.EasyCreate(words[len(words)-1])
			bamOut[i] = sam.NewBamWriter(bamOutHandle[i], headers[i])
		}
	}

	enclosingReads := make([][]*sam.Sam, len(inputFiles)) // first index is sample
	observedLengths := make([][]int, len(inputFiles))     // first index is sample
	var currVcf vcf.Vcf
	var i int
	alignerInput := make(chan sam.Sam, 100)
	alignerOutput := make(chan sam.Sam, 100)
	for j := 0; j < alignerThreads; j++ {
		ref = fasta.NewSeeker(refFile, "")
		defer ref.Close()
		go realign.RealignIndels(alignerInput, alignerOutput, ref)
	}
	mm := new(gmm.MixtureModel)
	gaussians := make([][]float64, 2)
	var floatSlice []float64
	var converged bool
	for _, region := range targets {
		for i = range inputFiles {
			enclosingReads[i], observedLengths[i] = getLenghtDist(enclosingReads[i], targetPadding, minMapQ, minFlankOverlap, removeDups, bamIdxs[i], region, br[i], bamOut[i], alignerInput, alignerOutput)
			if bamOutPfx != "" {
				for j := range enclosingReads[i] {
					sam.WriteToBamFileHandle(bamOut[i], *enclosingReads[i][j], 0)
				}
			}
			slices.Sort(observedLengths[i])
		}

		// TODO add idx for bulk sample
		converged = runMixtureModel(observedLengths[0], mm, &floatSlice)
		if !converged {
			continue
		}

		if debug > 0 {
			fmt.Println(region, len(observedLengths), printLengths(observedLengths))
			fmt.Println(mm.Means)
		}
		if debug > 1 {
			gaussians[0] = gaussianHist(mm.Weights[0], mm.Means[0], mm.Stdev[0])
			gaussians[1] = gaussianHist(mm.Weights[1], mm.Means[1], mm.Stdev[1])
			plot(observedLengths, minReads, gaussians)
		}
		//currVcf = callGenotypes(ref, region, minReads, enclosingReads, observedLengths, mm)
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
	err = vcfOut.Close()
	exception.PanicOnErr(err)
}

func callGenotypes(ref *fasta.Seeker, region bed.Bed, minReads int, enclosingReads [][]*sam.Sam, observedLengths [][]int, mm *gmm.MixtureModel) vcf.Vcf {
	var ans vcf.Vcf
	repeatUnitLen, refNumRepeats := parseRepeatSeq(region.Name)
	refRepeatLen := refNumRepeats * len(repeatUnitLen)
	ans.Chr = region.Chrom
	ans.Pos = region.ChromStart
	refSeq, err := fasta.SeekByName(ref, region.Chrom, region.ChromStart, region.ChromEnd)
	exception.PanicOnErr(err)
	dna.AllToUpper(refSeq)
	ans.Ref = dna.BasesToString(refSeq)
	//if len(ans.Ref) != refRepeatLen {
	//	log.Panicf("ERROR: %s ref seq is \n%s\n the length of %d does not match expected %d from bed file.", region, ans.Ref[1:], len(ans.Ref), refRepeatLen)
	//}

	ans.Id = region.Name
	altLens := make([]int, 2)
	var refLenDiff int
	for i, l := range mm.Means {
		altLens[i] = int(math.Round(l))
		refLenDiff = refRepeatLen - altLens[i]
		for _, alts := range ans.Alt {
			if len(alts) == altLens[i] {
				refLenDiff = 0 // to engage break below
			}
		}
		if refLenDiff == 0 {
			continue
		}
		ans.Alt = append(ans.Alt, ans.Ref[0:len(ans.Ref)-refLenDiff-1])
	}

	info := new(strings.Builder)
	info.WriteString(fmt.Sprintf("RefLength=%d", refRepeatLen))
	info.WriteString(";Means=")
	for i, j := range mm.Means {
		if i > 0 {
			info.WriteByte(',')
		}
		info.WriteString(fmt.Sprintf("%f", j))
	}
	info.WriteString(";Stdev=")
	for i, j := range mm.Stdev {
		if i > 0 {
			info.WriteByte(',')
		}
		info.WriteString(fmt.Sprintf("%f", j))
	}
	info.WriteString(";Weights=")
	for i, j := range mm.Weights {
		if i > 0 {
			info.WriteByte(',')
		}
		info.WriteString(fmt.Sprintf("%f", j))
	}
	ans.Info = info.String()

	fmt.Println(ans)
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
	realignReads(reads, minMapQ, alignerInput, alignerOutput) // read order in slice may change

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
	repeatSeq, _ := parseRepeatSeq(region.Name)
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
	return maxLength // TODO divide by repeat unit length???
}

func parseRepeatSeq(s string) ([]dna.Base, int) {
	var words []string
	if strings.Contains(s, "x") {
		words = strings.Split(s, "x")
	}
	num, err := strconv.Atoi(words[0])
	exception.PanicOnErr(err)
	return dna.StringToBases(words[1]), num
}

func dedup(reads []*sam.Sam) []*sam.Sam {
	for i := 1; i < len(reads); i++ {
		if reads[i].GetChromStart() == reads[i-1].GetChromStart() && reads[i].GetChromEnd() == reads[i-1].GetChromEnd() {
			slices.Delete(reads, i, i+1)
		}
	}
	return reads
}

// read order may change
func realignReads(reads []sam.Sam, minMapQ int, alignerInput chan<- sam.Sam, alignerOutput <-chan sam.Sam) {
	var readsSkipped, readsReceived int

	// count how many reads will be skipped over for realignment
	for i := range reads {
		if minMapQ != -1 && reads[i].MapQ < uint8(minMapQ) {
			readsSkipped++
		}
	}

	// start streaming reads to aligner
	go sendReads(reads, minMapQ, alignerInput)

	// start receiving aligned reads
	for read := range alignerOutput {
		reads[readsReceived] = read
		readsReceived++

		// break when all reads sent for alignment have been received
		if readsReceived+readsSkipped == len(reads) {
			reads = reads[0:readsReceived]
			break
		}
	}
}

func sendReads(reads []sam.Sam, minMapQ int, alignerInput chan<- sam.Sam) {
	for i := range reads {
		if minMapQ != -1 && reads[i].MapQ < uint8(minMapQ) {
			continue
		}
		alignerInput <- reads[i]
	}
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

func plot(observedLengths [][]int, minReads int, gaussians [][]float64) {
	readsPerSample := make([]int, len(observedLengths))
	p := make([][]float64, len(observedLengths))
	for i := range observedLengths {
		p[i] = make([]float64, 100)
		for j := range observedLengths[i] {
			p[i][observedLengths[i][j]]++
			readsPerSample[i]++
		}
	}
	if len(observedLengths) == 1 && readsPerSample[0] < minReads {
		return
	}
	fmt.Println(asciigraph.PlotMany(gaussians, asciigraph.Precision(0), asciigraph.SeriesColors(
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

	for i := range p {
		if readsPerSample[i] < minReads {
			continue
		}
		//if i != 0 {
		//	continue
		//}
		fmt.Println(asciigraph.Plot(p[i], asciigraph.Height(5), asciigraph.Precision(0), asciigraph.SeriesColors(asciigraph.AnsiColor(i))))
	}

	//fmt.Println(asciigraph.PlotMany(p, asciigraph.Precision(0), asciigraph.SeriesColors(
	//	asciigraph.Red,
	//	asciigraph.Yellow,
	//	asciigraph.Green,
	//	asciigraph.Blue,
	//	asciigraph.Cyan,
	//	asciigraph.BlueViolet,
	//	asciigraph.Brown,
	//	asciigraph.Gray,
	//	asciigraph.Orange,
	//	asciigraph.Olive,
	//), asciigraph.Height(10)))
}

func gaussianHist(weight, mean, stdev float64) []float64 {
	y := make([]float64, 100)
	for x := range y {
		y[x] = gaussianY(float64(x), weight, mean, stdev)
	}
	return y
}

func gaussianY(x, a, b, c float64) float64 {
	top := math.Pow(x-b, 2)
	bot := 2 * c * c
	return a * math.Exp(-top/bot)
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

func runMixtureModel(data []int, mm *gmm.MixtureModel, f *[]float64) bool {
	if cap(*f) >= len(data) {
		*f = (*f)[0:len(data)]
	} else {
		*f = make([]float64, len(data))
	}

	for i := range data {
		(*f)[i] = float64(data[i])
	}
	converged, _ := gmm.RunMixtureModel(*f, 2, 50, 50, mm)
	return converged
}
