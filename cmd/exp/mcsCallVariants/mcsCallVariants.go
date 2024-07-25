package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/call"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
)

func usage() {
	fmt.Print(
		"mcsCallVariants - Call variants from META-CS data processed with annotateReadFamilies.\n" +
			"Usage:\n" +
			"mcsCallVariants [options] -i input.bam -b input.bed -r reference.fasta > output.vcf\n\n")
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
	var excludeBeds inputFiles
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile")
	memprofile := flag.String("memprofile", "", "write memory profile")
	input := flag.String("i", "", "Input bam file. Must be indexed.")
	output := flag.String("o", "stdout", "Output VCF file.")
	bedFile := flag.String("b", "", "Input bed file with coordinates of read families, read family ID, and read counts for watson and crick strands. Generated with -bed option in annotateReadFamilies.")
	flag.Var(&excludeBeds, "e", "Bed file(s) with regions to exclude from analysis. May be declared more than once with additional -e flags. Strongly recommended to mask regions with poor mappability. Note that any family OVERLAPPING an excluded region will be removed from analysis.")
	ref := flag.String("r", "", "Fasta file with reference genome used to align input bam. Must be indexed.")
	totalDepth := flag.Int("a", 8, "Minimum total depth of read family for variant consideration.")
	strandedDepth := flag.Int("s", 4, "Minimum depth of independent watson and crick strands for variant consideration. When set to 0, caller runs in unstranded mode merging read counts from watson and crick strands.")
	endPad := flag.Int("ignoreEnds", 3, "Ignore bases within # of end of a read.")
	minMapQ := flag.Int("minMapQ", 20, "Minimum mapping quality.")
	minReadFamilyLength := flag.Int("minReadFamilyLength", 100, "Minimum length in bp of read family for inclusion in analysis. Empirical evidence suggests errors are more common in small fragments.")
	maxSoftClipFraction := flag.Float64("maxSoftClipFraction", 0.2, "Maximum fraction of read that may be soft clipped.")
	countOverlappingPairs := flag.Bool("countOverlappingPairs", false, "Count both reads in overlapping regions of read pairs. By only 1 base is contributed in overlapping regions of read pairs.")
	allowSuppAln := flag.Bool("allowSupplementaryAlignments", false, "Allow variants using reads that have supplementary alignments annotated.")
	minAf := flag.Float64("minAF", 0.9, "Minimum fraction of reads with alternate allele **Within a read family and within strand** to be considered a variant.")
	minBaseQuality := flag.Int("minBaseQuality", 30, "Minimum base quality to be considered for calling. Bases below threshold will be ignored.")
	baseQualPenalty := flag.Float64("baseQualPenalty", 0.5, "Penalty for positions with low quality base. Each read with a base < minBaseQuality counts towards baseQualPenalty fraction of a read for allele frequency calculations. Note that low quality bases are N-masked and so will always count AGAINST the alternate allele. (e.g. by default each read with a low quality base counts as 0.5 reads for allele frequency determination.")
	maxOverlappingFamilies := flag.Int("maxOverlappingFamilies", 20, "Maximum number of overlapping read families for site to be considered for calling. Low number avoids regions with many misalignments (e.g. centromeres) reducing memory usage. Set to -1 for no limit. Analyzed bed will be bedfile.analysis.bed")
	callSingleStrand := flag.Bool("ss", false, "Include single-stranded variants in output VCF. Single-stranded calling uses the same a and s minimum values as double-stranded calling but requires perfect asymmetry between strands such that 100% of reads carry the variant on strand 1 and 0% of reads carry the variant on strand 2. Single-stranded calls will have 'SS' in the INFO field.")
	minContigSize := flag.Int("minContigSize", 10_000_000, "Remove families mapping to contigs of length < minContigSize. The default value cuts out common decoy sequences and chrM from the human genome while keeping chr1-22,X,Y.")
	maxVariantsPerReadFamily := flag.Int("maxVariantsPerReadFamily", 3, "Maximum number of variants that are allowed to be called within a single read family. If a read family has more variants than this limit, all variants from the read family will be discarded.")
	threads := flag.Int("threads", 1, "Number of processor threads to use for calling. Output VCF will be out of order with threads > 1.")
	debugLevel := flag.Int("verbose", 0, "Level of verbosity in log.")
	debugOut := flag.String("debugLog", "", "Print debug logs to file. File may be large. Must be run with threads == 1 for coherent output. ")
	flag.Parse()

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

	if len(excludeBeds) == 0 {
		log.Println("WARNING: -e was not declared. It is strongly recommended to mask regions with poor mappability.")
	}

	if *threads == 0 {
		log.Fatal("ERROR: threads must be >= 1.")
	}

	if *input == "" || *bedFile == "" || *ref == "" {
		usage()
		log.Fatal("ERROR: must specify bam (-i), bed (-b), and fasta (-r).")
	}

	if *strandedDepth*2 > *totalDepth {
		log.Fatal("ERROR: -s * 2 should not be larger than -a")
	}

	call.Call(*input, *output, *ref, *bedFile, excludeBeds, uint8(*minMapQ), *totalDepth, *strandedDepth, *allowSuppAln, *minAf, *minBaseQuality, *minContigSize, *minReadFamilyLength, *baseQualPenalty, *maxSoftClipFraction, *endPad, *maxOverlappingFamilies, *countOverlappingPairs, *callSingleStrand, *maxVariantsPerReadFamily, *debugLevel, *threads, *debugOut)

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
