package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/call"
	"github.com/vertgenlab/gonomics/exception"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
)

func callUsage(callFlags *flag.FlagSet) {
	fmt.Print(
		"call - call variants from duplex data processed with 'duplextools pair'\n\n" +
			"Usage:\n" +
			"  duplextools call [options] -i input.bam -b input.bed -r reference.fasta > output.vcf\n\n" +
			"Options:\n")
	callFlags.PrintDefaults()
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

func runCall(args []string) {
	var err error
	callFlags := flag.NewFlagSet("call", flag.ExitOnError)

	var excludeBeds inputFiles
	cpuprofile := callFlags.String("cpuprofile", "", "write cpu profile")
	memprofile := callFlags.String("memprofile", "", "write memory profile")
	input := callFlags.String("i", "", "Input bam file. Must be indexed.")
	output := callFlags.String("o", "stdout", "Output VCF file.")
	bedFile := callFlags.String("b", "", "Input bed file with coordinates of read families, read family ID, and read counts for watson and crick strands. Generated with -bed option in annotateReadFamilies.")
	callFlags.Var(&excludeBeds, "e", "Bed file(s) with regions to exclude from analysis. May be declared more than once with additional -e flags. Strongly recommended to mask regions with poor mappability. Note that any family OVERLAPPING an excluded region will be removed from analysis.")
	ref := callFlags.String("r", "", "Fasta file with reference genome used to align input bam. Must be indexed.")
	totalDepth := callFlags.Int("a", 8, "Minimum total depth of read family for variant consideration.")
	strandedDepth := callFlags.Int("s", 4, "Minimum depth of independent watson and crick strands for variant consideration. When set to 0, caller runs in unstranded mode merging read counts from watson and crick strands.")
	endPad := callFlags.Int("ignoreEnds", 3, "Ignore bases within # of end of a read.")
	minMapQ := callFlags.Int("minMapQ", 20, "Minimum mapping quality.")
	minReadFamilyLength := callFlags.Int("minReadFamilyLength", 100, "Minimum length in bp of read family for inclusion in analysis. Empirical evidence suggests errors are more common in small fragments.")
	maxSoftClipFraction := callFlags.Float64("maxSoftClipFraction", 0.2, "Maximum fraction of read that may be soft clipped.")
	countOverlappingPairs := callFlags.Bool("countOverlappingPairs", false, "Count both reads in overlapping regions of read pairs. By only 1 base is contributed in overlapping regions of read pairs.")
	allowSuppAln := callFlags.Bool("allowSupplementaryAlignments", false, "Allow variants using reads that have supplementary alignments annotated.")
	minAf := callFlags.Float64("minAF", 0.9, "Minimum fraction of reads with alternate allele **Within a read family and within strand** to be considered a variant.")
	minBaseQuality := callFlags.Int("minBaseQuality", 30, "Minimum base quality to be considered for calling. Bases below threshold will be ignored.")
	baseQualPenalty := callFlags.Float64("baseQualPenalty", 0.5, "Penalty for positions with low quality base. Each read with a base < minBaseQuality counts towards baseQualPenalty fraction of a read for allele frequency calculations. Note that low quality bases are N-masked and so will always count AGAINST the alternate allele. (e.g. by default each read with a low quality base counts as 0.5 reads for allele frequency determination.")
	maxOverlappingFamilies := callFlags.Int("maxOverlappingFamilies", 20, "Maximum number of overlapping read families for site to be considered for calling. Low number avoids regions with many misalignments (e.g. centromeres) reducing memory usage. Set to -1 for no limit. Analyzed bed will be bedfile.analysis.bed")
	callSingleStrand := callFlags.Bool("ss", false, "Include single-stranded variants in output VCF. Single-stranded calling uses the same a and s minimum values as double-stranded calling but requires perfect asymmetry between strands such that 100% of reads carry the variant on strand 1 and 0% of reads carry the variant on strand 2. Single-stranded calls will have 'SS' in the INFO field.")
	minContigSize := callFlags.Int("minContigSize", 10_000_000, "Remove families mapping to contigs of length < minContigSize. The default value cuts out common decoy sequences and chrM from the human genome while keeping chr1-22,X,Y.")
	maxVariantsPerReadFamily := callFlags.Int("maxVariantsPerReadFamily", 3, "Maximum number of variants that are allowed to be called within a single read family. If a read family has more variants than this limit, all variants from the read family will be discarded.")
	threads := callFlags.Int("threads", 1, "Number of processor threads to use for calling. Output VCF will be out of order with threads > 1.")
	debugLevel := callFlags.Int("verbose", 0, "Level of verbosity in log.")
	debugOut := callFlags.String("debugLog", "", "Print debug logs to file. File may be large. Must be run with threads == 1 for coherent output. ")

	err = callFlags.Parse(args)
	exception.PanicOnErr(err)
	callFlags.Usage = func() { callUsage(callFlags) }

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			errExit(err.Error())
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			errExit(err.Error())
		}
		defer pprof.StopCPUProfile()
	}

	if len(excludeBeds) == 0 {
		log.Println("WARNING: -e was not declared. It is strongly recommended to mask regions with poor mapability.")
	}

	if *threads == 0 {
		callFlags.Usage()
		errExit("\nERROR: threads must be >= 1")
	}

	if *input == "" || *bedFile == "" || *ref == "" {
		callFlags.Usage()
		errExit("\nERROR: must specify bam (-i), bed (-b), and fasta (-r)")
	}

	if *strandedDepth*2 > *totalDepth {
		callFlags.Usage()
		errExit("\nERROR: -s * 2 should not be larger than -a")
	}

	call.Call(*input, *output, *ref, *bedFile, excludeBeds, uint8(*minMapQ), *totalDepth, *strandedDepth, *allowSuppAln, *minAf, *minBaseQuality, *minContigSize, *minReadFamilyLength, *baseQualPenalty, *maxSoftClipFraction, *endPad, *maxOverlappingFamilies, *countOverlappingPairs, *callSingleStrand, *maxVariantsPerReadFamily, *debugLevel, *threads, *debugOut)

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			errExit(err.Error())
		}
		defer f.Close()
		runtime.GC() // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			errExit(err.Error())
		}
	}
}
