package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/extract"
	"github.com/vertgenlab/gonomics/exception"
)

func extractUsage(extractFlags *flag.FlagSet) {
	fmt.Print(
		"extract - Process raw FASTQ files to an unmapped BAM with barcode tags\n\n" +
			"Usage:\n" +
			"  duplextools extract [options] -1 r1.fq.gz -2 r2.fq.gz > output.bam\n\n" +
			"Options:\n")
	extractFlags.PrintDefaults()
}

func runExtract(args []string) {
	var err error
	extractFlags := flag.NewFlagSet("extract", flag.ExitOnError)

	r1 := extractFlags.String("1", "", "FASTQ file containing R1 reads. May be gzipped.")
	r2 := extractFlags.String("2", "", "FASTQ file containing R2 reads. May be gzipped.")
	outfile := extractFlags.String("o", "stdout", "Output BAM file.")
	missingBcFile := extractFlags.String("missing", "", "Output BAM file for records with missing barcodes")

	err = extractFlags.Parse(args)
	exception.PanicOnErr(err)
	extractFlags.Usage = func() { extractUsage(extractFlags) }

	if *r1 == "" || *r2 == "" {
		extractFlags.Usage()
		errExit("\nERROR: must have inputs for -1 and -2")
	}

	extract.Extract(*r1, *r2, *outfile, *missingBcFile)
}
