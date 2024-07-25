package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/extract"
	"log"
)

func usage() {
	fmt.Print(
		"mcsFqToBam - Process raw FASTQ files generated with META-CS to an unmapped BAM with barcode tags.\n\n" +
			"Usage:\n" +
			"  mcsFqToBam [options] -1 r1.fq.gz -2 r2.fq.gz -o output.bam\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func main() {
	r1 := flag.String("1", "", "FASTQ file containing R1 reads. May be gzipped.")
	r2 := flag.String("2", "", "FASTQ file containing R2 reads. May be gzipped.")
	outfile := flag.String("o", "stdout", "Output BAM file.")
	missingBcFile := flag.String("missing", "", "Output BAM file for records with missing barcodes")
	flag.Parse()
	flag.Usage = usage

	if *r1 == "" || *r2 == "" {
		flag.Usage()
		log.Fatalln("ERROR: must input read1 file and read2 file")
	}

	extract.Extract(*r1, *r2, *outfile, *missingBcFile)
}
