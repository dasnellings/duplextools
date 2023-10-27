package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/fai"
	"github.com/dasnellings/duplexTools/gmm"
	"github.com/dasnellings/duplexTools/realign"
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
	"gonum.org/v1/gonum/stat"
	"io"
	"log"
	"math"
	"os"
	"path"
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
	var inputDir *string = flag.String("inputDir", "", "Directory with BAM files to be used as inputs. Uses all files in the directory ending with \".bam\". Can be used instead of -i.")
	var ref *string = flag.String("r", "", "Reference genome. Must be the same reference used for generating the BAM file.")
	var targets *string = flag.String("t", "", "BED file of targeted repeats. The 4th column must be the sequence of one repeat unit (e.g. CA for a CACACACA repeat), or 'RepeatLen'x'RepeatSeq' (e.g. 10xCA).")
	var output *string = flag.String("o", "stdout", "Output VCF file.")
	var lenOut *string = flag.String("lenOut", "", "Output a bed file with additional columns for determined read lengths for each sample.")
	var bamOut *string = flag.String("bamOutPfx", "", "Output a BAM file with realigned reads. Only outputs reads that inform called genotypes. File will be named 'bamOutPfx'_'originalFilename'.")
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
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	if *inputDir != "" {
		inputs = getInputsFromDir(*inputDir)
	}

	if len(inputs) == 0 || *ref == "" {
		usage()
		log.Fatalln("ERROR: must input a BAM file with -i")
	}

	debug = *debugVal

	if *minMapQ > math.MaxUint8 {
		log.Fatalf("minMapQ out of range. max: %d\n", math.MaxUint8)
	}

	genotypeTargetRepeats(inputs, *ref, *targets, *output, *bamOut, *lenOut, *targetPadding, *minFlankOverlap, *minMapQ, *minReads, !*allowDups, *alignerThreads)

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

func getInputsFromDir(dir string) []string {
	var inputs []string
	files, err := os.ReadDir(dir)
	exception.PanicOnErr(err)
	for i := range files {
		if strings.HasSuffix(files[i].Name(), ".bam") {
			inputs = append(inputs, files[i].Name())
		}
	}
	return inputs
}

func genotypeTargetRepeats(inputFiles []string, refFile, targetsFile, outputFile, bamOutPfx, lenOutFile string, targetPadding, minFlankOverlap, minMapQ, minReads int, removeDups bool, alignerThreads int) {
	var err error
	var ref *fasta.Seeker
	var lenOut *fileio.EasyWriter
	buf := new([2][11]float64)
	readBuf := new([]float64)
	targets := bed.Read(targetsFile)
	vcfOut := fileio.EasyCreate(outputFile)
	defer cleanup(vcfOut)
	vcfHeader := generateVcfHeader(strings.Join(inputFiles, "\t"), refFile)
	vcf.NewWriteHeader(vcfOut, vcfHeader)

	// get bam reader for each file
	br := make([]*sam.BamReader, len(inputFiles))
	headers := make([]sam.Header, len(inputFiles))
	bamIdxs := make([]sam.Bai, len(inputFiles))
	for i := range inputFiles {
		br[i], headers[i] = sam.OpenBam(inputFiles[i])
		defer cleanup(br[i])
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
			defer cleanup(bamOutHandle[i])
			defer cleanup(bamOut[i])
		}
	}

	if lenOutFile != "" {
		lenOut = fileio.EasyCreate(lenOutFile)
		fmt.Fprintf(lenOut, "#CHROM\tSTART\tEND\tREPEAT\t%s\n", strings.Join(inputFiles, "\t"))
		defer cleanup(lenOut)
	}

	enclosingReads := make([][]*sam.Sam, len(inputFiles)) // first index is sample
	observedLengths := make([][]int, len(inputFiles))     // first index is sample
	var currVcf vcf.Vcf
	alignerInput := make(chan sam.Sam, 1000)
	alignerOutput := make(chan sam.Sam, 1000)
	for j := 0; j < alignerThreads; j++ {
		ref = fasta.NewSeeker(refFile, "")
		defer cleanup(ref)
		go realign.RealignIndels(alignerInput, alignerOutput, ref)
	}

	mm := make([]*gmm.MixtureModel, len(inputFiles))
	tmpMm := make([]*gmm.MixtureModel, len(inputFiles))
	for i := 0; i < len(inputFiles); i++ {
		mm[i] = new(gmm.MixtureModel)
		tmpMm[i] = new(gmm.MixtureModel)
	}

	gaussians := make([][]float64, 2)
	var floatSlices [][]float64 = make([][]float64, len(inputFiles))
	var converged, anyConverged, passingVariant bool
	var repeatUnit []dna.Base
	for _, region := range targets {
		repeatUnit, _ = parseRepeatSeq(region.Name)
		anyConverged = false
		for i := range inputFiles {
			enclosingReads[i], observedLengths[i] = getLenghtDist(enclosingReads[i], targetPadding, minMapQ, minFlankOverlap, removeDups, bamIdxs[i], region, br[i], bamOut[i], alignerInput, alignerOutput)
			if bamOutPfx != "" {
				for j := range enclosingReads[i] {
					sam.WriteToBamFileHandle(bamOut[i], *enclosingReads[i][j], 0)
				}
			}
			slices.Sort(observedLengths[i])

			converged, tmpMm[i], mm[i] = runMixtureModel(observedLengths[i], tmpMm[i], mm[i], &floatSlices[i])
			if converged {
				anyConverged = true
			}
		}

		if !anyConverged {
			continue
		}

		if lenOut != nil {
			fmt.Fprintf(lenOut, "%s%s\n", bed.ToString(region, 4), printLengths(observedLengths))
		}

		if debug > 0 {
			//val, counts := sliceToCounts(mm[0].Data)
			//for i := range val {
			//	fmt.Printf("%d:%d\t", int(val[i]), counts[i])
			//}
			//fmt.Println()
			for i := range mm {
				for k := range mm[i].Means {
					fmt.Printf("k=%d mu=%0.2f stdev=%0.2f\tloglikelihood=%0.4g\n", k, mm[i].Means[k], mm[i].Stdev[k], mm[i].LogLikelihood)
					testPulseFitKS(mm[i], k, len(repeatUnit), buf, readBuf, true)
					testPulseFitHeuristic(mm[i], k, len(repeatUnit), true)
				}
			}
			plot(observedLengths, minReads, mm, gaussians)
		}

		currVcf, passingVariant = callGenotypes(ref, region, minReads, enclosingReads, observedLengths, mm, buf, readBuf)
		if passingVariant {
			vcf.WriteVcf(vcfOut, currVcf)
		}
	}
	close(alignerInput)
	close(alignerOutput)
}

func callGenotypes(ref *fasta.Seeker, region bed.Bed, minReads int, enclosingReads [][]*sam.Sam, observedLengths [][]int, mm []*gmm.MixtureModel, buf *[2][11]float64, readBuf *[]float64) (vcf.Vcf, bool) {
	var ans vcf.Vcf
	repeatUnitLen, refNumRepeats := parseRepeatSeq(region.Name)
	refRepeatLen := refNumRepeats * len(repeatUnitLen)
	ans.Chr = region.Chrom
	ans.Pos = region.ChromStart
	refSeq, err := fasta.SeekByName(ref, region.Chrom, region.ChromStart, region.ChromEnd)
	exception.PanicOnErr(err)
	dna.AllToUpper(refSeq)
	ans.Ref = dna.BasesToString(refSeq)
	ans.Ref = "*" // TODO Remove
	//if len(ans.Ref) != refRepeatLen {
	//	log.Panicf("ERROR: %s ref seq is \n%s\n the length of %d does not match expected %d from bed file.", region, ans.Ref[1:], len(ans.Ref), refRepeatLen)
	//}

	ans.Id = region.Name

	/*
		altLens := make([]int, 2)
		var refLenDiff int
		for i, l := range mm[0].Means {
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
	*/
	ans.Alt = append(ans.Alt, "*")
	ans.Filter = "."
	ans.Id = region.Name
	ans.Format = []string{"GT", "DP", "MU", "SD", "WT", "LL", "AD", "KS", "CG", "HS", "HG", "RL"}
	ans.Samples = make([]vcf.Sample, len(mm))
	var goodnessOfFit0, goodnessOfFit1, pulseHeuristic0, pulseHeuristic1 float64
	var allele0Reads, allele1Reads, minKsLen0, minKsLen1, optimalHeuristicLen0, optimalHeuristicLen1 int
	var readLenString0, readLenString1 string

	//for j := range mm[0].Data {
	//	fmt.Printf("%0.0f, %0.1f, %0.1f\t", mm[0].Data[j], mm[0].Posteriors[0][j], mm[0].Posteriors[1][j])
	//}

	for i := range ans.Samples {
		ans.Samples[i].FormatData = make([]string, 12)
		ans.Samples[i].FormatData[1] = fmt.Sprintf("%d", len(observedLengths[i]))
		if mm[i].LogLikelihood == math.MaxFloat64 {
			ans.Samples[i].FormatData[2] = "."
			ans.Samples[i].FormatData[3] = "."
			ans.Samples[i].FormatData[4] = "."
			ans.Samples[i].FormatData[5] = "."
			ans.Samples[i].FormatData[6] = "."
			ans.Samples[i].FormatData[7] = "."
			ans.Samples[i].FormatData[8] = "."
			ans.Samples[i].FormatData[9] = "."
			ans.Samples[i].FormatData[10] = "."
			ans.Samples[i].FormatData[11] = "."
			continue
		}
		ans.Samples[i].FormatData[5] = fmt.Sprintf("%.1g", mm[i].LogLikelihood)

		goodnessOfFit0, allele0Reads, minKsLen0 = testPulseFitKS(mm[i], 0, len(repeatUnitLen), buf, readBuf, false)
		goodnessOfFit1, allele1Reads, minKsLen1 = testPulseFitKS(mm[i], 1, len(repeatUnitLen), buf, readBuf, false)
		pulseHeuristic0, _, optimalHeuristicLen0 = testPulseFitHeuristic(mm[i], 0, len(repeatUnitLen), false)
		pulseHeuristic1, _, optimalHeuristicLen1 = testPulseFitHeuristic(mm[i], 1, len(repeatUnitLen), false)
		readLenString0 = getRunLengthEncoding(getReadsForK(mm[i], 0, readBuf))
		readLenString1 = getRunLengthEncoding(getReadsForK(mm[i], 1, readBuf))

		if mm[i].Means[0] < mm[i].Means[1] {
			ans.Samples[i].FormatData[2] = fmt.Sprintf("%.1f,%.1f", mm[i].Means[0], mm[i].Means[1])
			ans.Samples[i].FormatData[3] = fmt.Sprintf("%.1f,%.1f", mm[i].Stdev[0], mm[i].Stdev[1])
			ans.Samples[i].FormatData[4] = fmt.Sprintf("%.1f,%.1f", mm[i].Weights[0], mm[i].Weights[1])
			ans.Samples[i].FormatData[6] = fmt.Sprintf("%d,%d", allele0Reads, allele1Reads)
			ans.Samples[i].FormatData[7] = fmt.Sprintf("%.3f,%.3f", goodnessOfFit0, goodnessOfFit1)
			ans.Samples[i].FormatData[8] = fmt.Sprintf("%d,%d", minKsLen0, minKsLen1)
			ans.Samples[i].FormatData[9] = fmt.Sprintf("%.3f,%.3f", pulseHeuristic0, pulseHeuristic1)
			ans.Samples[i].FormatData[10] = fmt.Sprintf("%d,%d", optimalHeuristicLen0, optimalHeuristicLen1)
			ans.Samples[i].FormatData[11] = fmt.Sprintf("%s;%s", readLenString0, readLenString1)
		} else {
			ans.Samples[i].FormatData[2] = fmt.Sprintf("%.1f,%.1f", mm[i].Means[1], mm[i].Means[0])
			ans.Samples[i].FormatData[3] = fmt.Sprintf("%.1f,%.1f", mm[i].Stdev[1], mm[i].Stdev[0])
			ans.Samples[i].FormatData[4] = fmt.Sprintf("%.1f,%.1f", mm[i].Weights[1], mm[i].Weights[0])
			ans.Samples[i].FormatData[6] = fmt.Sprintf("%d,%d", allele1Reads, allele0Reads)
			ans.Samples[i].FormatData[7] = fmt.Sprintf("%.3f,%.3f", goodnessOfFit1, goodnessOfFit0)
			ans.Samples[i].FormatData[8] = fmt.Sprintf("%d,%d", minKsLen1, minKsLen0)
			ans.Samples[i].FormatData[9] = fmt.Sprintf("%.3f,%.3f", pulseHeuristic1, pulseHeuristic0)
			ans.Samples[i].FormatData[10] = fmt.Sprintf("%d,%d", optimalHeuristicLen1, optimalHeuristicLen0)
			ans.Samples[i].FormatData[11] = fmt.Sprintf("%s;%s", readLenString1, readLenString0)
		}
	}

	ans.Info = fmt.Sprintf("RefLength=%d", refRepeatLen)
	return ans, true
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

func generateVcfHeader(samples string, referenceFile string) vcf.Header {
	var header vcf.Header
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, fmt.Sprintf("##reference=%s", path.Clean(referenceFile)))
	header.Text = append(header.Text, strings.TrimSuffix(fai.IndexToVcfHeader(fai.ReadIndex(referenceFile+".fai")), "\n"))
	header.Text = append(header.Text, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	header.Text = append(header.Text, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=MU,Number=2,Type=Float,Description=\"Mean repeat length of each allele determined by gaussian mixture modelling.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=SD,Number=2,Type=Float,Description=\"Standard deviation of the repeat length of each allele determined by gaussian mixture modelling.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=WT,Number=2,Type=Float,Description=\"Weight assigned to each allele (rough estimate of allele frequency) determined by gaussian mixture modelling.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=LL,Number=1,Type=Float,Description=\"Negative log likelihood of gaussian mixture model.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Number of reads assigned to each allele based on posteriors from gaussian modelling.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=KS,Number=2,Type=Float,Description=\"Kolmogorov-Smirnov (KS) statistic for fit of data to oscillating slippage model dependent on repeat unit length.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=CG,Number=2,Type=Integer,Description=\"Optimal repeat length fit as determined by minimum KS statistic.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=HS,Number=2,Type=Float,Description=\"Heuristic score for fit of data to oscillating slippage model dependent on repeat unit length. Higher values indicate better fit to slippage model\">")
	header.Text = append(header.Text, "##FORMAT=<ID=HG,Number=2,Type=Integer,Description=\"Optimal repeat length fit as determined by maximum heuristic score.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=RL,Number=2,Type=String,Description=\"Run length encoding of read lengths for each allele separated by semicolons.\">")
	header.Text = append(header.Text, "##INFO=<ID=RefLength,Number=1,Type=Integer,Description=\"Length in bp of the repeat in the reference genome.\">")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", strings.Replace(samples, ".bam", "", -1)))
	return header
}

func sliceToCounts(s []float64) (val []float64, count []int) {
	m := make(map[float64]int)
	for i := range s {
		m[s[i]]++
	}
	for k := range m {
		val = append(val, k)
	}

	slices.Sort(val)
	for i := range val {
		count = append(count, m[val[i]])
	}
	return
}

func plot(observedLengths [][]int, minReads int, mm []*gmm.MixtureModel, gaussians [][]float64) {
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

	for i := range p {
		if readsPerSample[i] < minReads {
			continue
		}
		//if i != 0 {
		//	continue
		//}
		fmt.Println(asciigraph.Plot(p[i], asciigraph.Height(5), asciigraph.Precision(0), asciigraph.SeriesColors(asciigraph.AnsiColor(i))))

		gaussians[0] = gaussianHist(mm[i].Weights[0], mm[i].Means[0], mm[i].Stdev[0])
		gaussians[1] = gaussianHist(mm[i].Weights[1], mm[i].Means[1], mm[i].Stdev[1])

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

func gaussianY(x, weight, mean, stdev float64) float64 {
	top := math.Pow(x-mean, 2)
	bot := 2 * stdev * stdev
	return weight * math.Exp(-top/bot)
}

func printLengths(a [][]int) string {
	if len(a) == 0 {
		return ""
	}
	s := new(strings.Builder)
	for i := range a {
		if len(a[i]) == 0 {
			s.WriteString("\tNA")
			continue
		}
		s.WriteString(fmt.Sprintf("\t%d", a[i][0]))
		for j := 1; j < len(a[i]); j++ {
			s.WriteString(fmt.Sprintf(",%d", a[i][j]))
		}
	}
	return s.String()
}

func runMixtureModel(data []int, mm, bestMm *gmm.MixtureModel, f *[]float64) (converged bool, newMm, newBestMm *gmm.MixtureModel) {
	if cap(*f) >= len(data) {
		*f = (*f)[0:len(data)]
	} else {
		*f = make([]float64, len(data))
	}

	for i := range data {
		(*f)[i] = float64(data[i])
	}

	for i := 0; i < 10; i++ {
		converged, _ = gmm.RunMixtureModel(*f, 2, 50, 50, mm)
		if i == 0 {
			mm, bestMm = bestMm, mm
			continue
		}
		if mm.LogLikelihood < bestMm.LogLikelihood {
			mm, bestMm = bestMm, mm
		}
	}
	return converged, mm, bestMm
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}

func testPulseFitHeuristic(mm *gmm.MixtureModel, k, period int, print bool) (float64, int, int) {
	var readsIncluded, readLen int
	var maxFitSum, startFitSum, lessFitSum, moreFitSum, less2FitSum, more2FitSum, startVal, lessVal, moreVal, less2Val, more2Val float64

	// test original peak as well as +/- 1 to account for slip-up/down bias
	var startPeak, lessPeak, morePeak, bestPeak, less2Peak, more2Peak int
	startPeak = int(math.Round(mm.Means[k]))
	lessPeak = startPeak - 1
	less2Peak = startPeak - 2
	morePeak = startPeak + 1
	more2Peak = startPeak + 2

	for i := range mm.Data {
		readLen = int(mm.Data[i])
		if getMaxK(mm.Posteriors, i) != k {
			continue
		}
		readsIncluded++
		startVal = gaussianY(mm.Data[i], 1, float64(startPeak), mm.Stdev[k]) // TODO values could be cached
		lessVal = gaussianY(mm.Data[i], 1, float64(lessPeak), mm.Stdev[k])
		moreVal = gaussianY(mm.Data[i], 1, float64(morePeak), mm.Stdev[k])
		less2Val = gaussianY(mm.Data[i], 1, float64(less2Peak), mm.Stdev[k])
		more2Val = gaussianY(mm.Data[i], 1, float64(more2Peak), mm.Stdev[k])

		if (startPeak-readLen)%period == 0 {
			startFitSum += startVal
		} else {
			startFitSum -= 1
		}

		if (lessPeak-readLen)%period == 0 {
			lessFitSum += lessVal
		} else {
			lessFitSum -= 1
		}

		if (morePeak-readLen)%period == 0 {
			moreFitSum += moreVal
		} else {
			moreFitSum -= 1
		}

		if (less2Peak-readLen)%period == 0 {
			less2FitSum += less2Val
		} else {
			less2FitSum -= 1
		}

		if (more2Peak-readLen)%period == 0 {
			more2FitSum += more2Val
		} else {
			more2FitSum -= 1
		}
	}

	maxFitSum = math.Max(startFitSum, math.Max(lessFitSum, moreFitSum))
	switch maxFitSum {
	case startFitSum:
		bestPeak = startPeak
	case lessFitSum:
		bestPeak = lessPeak
	case moreFitSum:
		bestPeak = morePeak
	}

	if print {
		log.Printf("Pulse fit heuristic:\tlen=%d\tvalue=%0.2f\treads=%d\n", lessPeak, ((lessFitSum/float64(readsIncluded))+1)/2, readsIncluded)
		log.Printf("Pulse fit heuristic:\tlen=%d\tvalue=%0.2f\treads=%d\n", startPeak, ((startFitSum/float64(readsIncluded))+1)/2, readsIncluded)
		log.Printf("Pulse fit heuristic:\tlen=%d\tvalue=%0.2f\treads=%d\n", morePeak, ((moreFitSum/float64(readsIncluded))+1)/2, readsIncluded)
	}
	return ((maxFitSum / float64(readsIncluded)) + 1) / 2, readsIncluded, bestPeak
}

func testPulseFitKS(mm *gmm.MixtureModel, k, period int, buf *[2][11]float64, readBuf *[]float64, print bool) (float64, int, int) {
	// test original peak as well as +/- 1 to account for slip-up/down bias
	var startPeak, lessPeak, morePeak, less2Peak, more2Peak int
	startPeak = int(math.Round(mm.Means[k]))
	lessPeak = startPeak - 1
	less2Peak = startPeak - 2
	morePeak = startPeak + 1
	more2Peak = startPeak + 2

	reads := getReadsForK(mm, k, readBuf)
	slices.Sort(reads)

	expectedK0less2Val, expectedK0less2Weight := getExpectedValuesForK(mm, k, buf, less2Peak, period, true)
	ansLess2 := stat.KolmogorovSmirnov(reads, nil, expectedK0less2Val, expectedK0less2Weight)
	if print {
		//fmt.Println(reads)
		//fmt.Println(expectedK0lessVal)
		//fmt.Println(expectedK0lessWeight)
		fmt.Printf("k=%d reads=%d peak=%d, ks=%0.4f\n", k, len(reads), less2Peak, ansLess2)
	}

	expectedK0lessVal, expectedK0lessWeight := getExpectedValuesForK(mm, k, buf, lessPeak, period, true)
	ansLess := stat.KolmogorovSmirnov(reads, nil, expectedK0lessVal, expectedK0lessWeight)
	if print {
		//fmt.Println(reads)
		//fmt.Println(expectedK0lessVal)
		//fmt.Println(expectedK0lessWeight)
		fmt.Printf("k=%d reads=%d peak=%d, ks=%0.4f\n", k, len(reads), lessPeak, ansLess)
	}

	expectedK0startVal, expectedK0startWeight := getExpectedValuesForK(mm, k, buf, startPeak, period, false)
	ansStart := stat.KolmogorovSmirnov(reads, nil, expectedK0startVal, expectedK0startWeight)
	if print {
		//fmt.Println(expectedK0startVal)
		//fmt.Println(expectedK0startWeight)
		fmt.Printf("k=%d reads=%d peak=%d, ks=%0.4f\n", k, len(reads), startPeak, ansStart)
	}

	expectedK0moreVal, expectedK0moreWeight := getExpectedValuesForK(mm, k, buf, morePeak, period, false)
	ansMore := stat.KolmogorovSmirnov(reads, nil, expectedK0moreVal, expectedK0moreWeight)
	if print {
		//fmt.Println(expectedK0moreVal)
		//fmt.Println(expectedK0moreWeight)
		fmt.Printf("k=%d reads=%d peak=%d, ks=%0.4f\n", k, len(reads), morePeak, ansMore)
	}

	expectedK0more2Val, expectedK0more2Weight := getExpectedValuesForK(mm, k, buf, more2Peak, period, false)
	ansMore2 := stat.KolmogorovSmirnov(reads, nil, expectedK0more2Val, expectedK0more2Weight)
	if print {
		//fmt.Println(expectedK0moreVal)
		//fmt.Println(expectedK0moreWeight)
		fmt.Printf("k=%d reads=%d peak=%d, ks=%0.4f\n", k, len(reads), more2Peak, ansMore2)
	}

	minScore := minslice(ansStart, ansLess, ansMore, ansMore2, ansLess2)
	switch minScore {
	case ansStart:
		return ansStart, len(reads), startPeak
	case ansLess:
		return ansLess, len(reads), lessPeak
	case ansLess2:
		return ansLess2, len(reads), less2Peak
	case ansMore:
		return ansMore, len(reads), morePeak
	case ansMore2:
		return ansMore2, len(reads), more2Peak
	default:
		panic("unreachable")
		return 6, 6, 6
	}
}

func getReadsForK(mm *gmm.MixtureModel, k int, readBuf *[]float64) []float64 {
	*readBuf = (*readBuf)[:0]
	if cap(*readBuf) < len(mm.Data) {
		*readBuf = make([]float64, 0, len(mm.Data))
	}

	for i := range mm.Data {
		if getMaxK(mm.Posteriors, i) != k {
			continue
		}
		*readBuf = append(*readBuf, mm.Data[i])
	}
	return *readBuf
}

func getExpectedValuesForK(mm *gmm.MixtureModel, k int, buf *[2][11]float64, peak, period int, updateWeights bool) ([]float64, []float64) {
	var currLen int = peak - (5 * period)
	for i := 0; i < 11; i++ {
		(*buf)[0][i] = float64(currLen)
		if updateWeights {
			(*buf)[1][i] = gaussianY(float64(currLen), 1, float64(peak), mm.Stdev[k])
		}
		currLen += period
	}
	return (*buf)[0][:], (*buf)[1][:]
}

func getMaxK(posteriors [][]float64, i int) int {
	var maxK int
	var maxVal float64
	for k := range posteriors {
		if posteriors[k][i] > maxVal {
			maxK = k
			maxVal = posteriors[k][i]
		}
	}
	return maxK
}

func min(a, b float64) float64 {
	if a < b {
		return a
	} else {
		return b
	}
}

func minslice(vals ...float64) float64 {
	minval := vals[0]
	for i := 1; i < len(vals); i++ {
		if vals[i] < minval {
			minval = vals[i]
		}
	}
	return minval
}

func getRunLengthEncoding(s []float64) string {
	ans := new(strings.Builder)
	var currVal float64
	var currCount int
	for i := range s {
		if s[i] != currVal {
			if ans.Len() == 0 && currVal != 0 {
				ans.WriteString(fmt.Sprintf("%d=%d", int(currVal), currCount))
			} else if currVal != 0 {
				ans.WriteString(fmt.Sprintf(",%d=%d", int(currVal), currCount))
			}
			currVal = s[i]
			currCount = 1
		} else {
			currCount++
		}
	}

	if ans.Len() == 0 && currVal != 0 {
		ans.WriteString(fmt.Sprintf("%d=%d", int(currVal), currCount))
	} else if currVal != 0 {
		ans.WriteString(fmt.Sprintf(",%d=%d", int(currVal), currCount))
	}
	return ans.String()
}
