package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/filter"
	"github.com/vertgenlab/gonomics/exception"
	"log"
)

func filterUsage(filterFlags *flag.FlagSet) {
	fmt.Print(
		"filter - filter constitutional (germline) variants from a vcf based on data from bulk sequencing\n\n" +
			"Usage:\n" +
			"  duplextools filter [options] -i input.vcf -b bulkGenomic.bam > output.vcf\n\n" +
			"Options:\n")
	filterFlags.PrintDefaults()
}

func runFilter(args []string) {
	var err error
	filterFlags := flag.NewFlagSet("filter", flag.ExitOnError)

	input := filterFlags.String("i", "", "Input VCF file with variant calls.")
	genomicVcf := filterFlags.String("g", "", "VCF file with germline variants from bulk sequencing.")
	genomicBam := filterFlags.String("b", "", "BAM file from bulk tissue. Must be indexed (.bai).")
	snpVcf := filterFlags.String("e", "", "VCF file with sites with known SNPs to exclude from analysis.")
	minCoverage := filterFlags.Int("minCoverage", 10, "Minimum coverage in bulk bam for consideration in output.")
	maxReadFrac := filterFlags.Float64("maxReadFrac", 0.1, "Maximum fraction of reads (minimum 1) in bulk sample for variant to be considered for output.")
	maxReads := filterFlags.Int("maxReads", 100000, "Maximum number of reads with alternate allele present in bulk sample to escape filtering (e.g. set to 1 to exclude all variants with >1 read with alternate allele in bulk sample")
	minBaseQuality := filterFlags.Int("minBaseQuality", 0, "Minimum base quality to be considered for calling. Bases below threshold will be ignored.")
	output := filterFlags.String("o", "stdout", "Output VCF file.")

	err = filterFlags.Parse(args)
	exception.PanicOnErr(err)
	filterFlags.Usage = func() { filterUsage(filterFlags) }

	if *input == "" || *genomicBam == "" {
		filterFlags.Usage()
		errExit("\nERROR: must have inputs for -i, and -b")
	}

	if *genomicVcf == "" {
		log.Println("WARNING: use of -g is STRONGLY RECOMMENDED if you are analyzing indels. It is useful, but not critical for analysing SNVs.")
	}

	filter.Filter(*input, *output, *genomicBam, *genomicVcf, *snpVcf, *minCoverage, *maxReadFrac, *maxReads, *minBaseQuality)
}
