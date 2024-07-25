package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/filter"
	"log"
)

func usage() {
	fmt.Print(
		"filterGermline - Filter germline variants from a vcf based on data from bulk sequencing.\n" +
			"Usage:\n" +
			"filterGermline [options] -i input.vcf -b bulkGenomic.bam > output.vcf\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input VCF file with variant calls.")
	genomicVcf := flag.String("g", "", "VCF file with germline variants from bulk sequencing.")
	genomicBam := flag.String("b", "", "BAM file from bulk tissue. Must be indexed (.bai).")
	snpVcf := flag.String("e", "", "VCF file with sites with known SNPs to exclude from analysis.")
	minCoverage := flag.Int("minCoverage", 10, "Minimum coverage in bulk bam for consideration in output.")
	maxReadFrac := flag.Float64("maxReadFrac", 0.1, "Maximum fraction of reads (minimum 1) in bulk sample for variant to be considered for output.")
	maxReads := flag.Int("maxReads", 100000, "Maximum number of reads with alternate allele present in bulk sample to escape filtering (e.g. set to 1 to exclude all variants with >1 read with alternate allele in bulk sample")
	minBaseQuality := flag.Int("minBaseQuality", 0, "Minimum base quality to be considered for calling. Bases below threshold will be ignored.")
	output := flag.String("o", "stdout", "Output VCF file.")
	flag.Parse()

	if *genomicVcf == "" {
		log.Println("WARNING: use of -g is STRONGLY RECOMMENDED if you are analyzing indels. It is useful, but not critical for analysing SNVs.")
	}

	if *input == "" || *genomicBam == "" {
		usage()
		log.Fatalln("ERROR: must have inputs for -i, and -b")
	}

	filter.Filter(*input, *output, *genomicBam, *genomicVcf, *snpVcf, *minCoverage, *maxReadFrac, *maxReads, *minBaseQuality)
}
