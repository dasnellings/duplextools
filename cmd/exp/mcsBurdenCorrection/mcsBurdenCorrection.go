package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/burden"
	"log"
)

func usage() {
	fmt.Print(
		"mcsBurdenCorrection - Correct META-CS mutation burden for trinucleotide (or other) context bias and allele count.\n" +
			"Applies method used for NanoSeq (for more information see Abascal et al. 2021 PMID: 33911282)\n" +
			"Usage:\n" +
			"mcsBurdenCorrection [options] -i input.vcf -b readFamilies.bed -r reference.fasta > summary.tsv\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input VCF file with final variant calls (somatic only, after filters).")
	ref := flag.String("r", "", "Reference FASTA file. Must be indexed (.fai).")
	bedfile := flag.String("b", "", "Bed file of read families used for calling. If ends were trimmed from read families for variant calling, it should be reflected in the input bed regions or bTrim. DO NOT INCLUDE FAMILIES THAT WERE EXCLUDED FROM ANALYSIS!!!")
	pad := flag.Int("pad", 1, "Number of context bases on either side of variant (e.g. 0 == T, 1 == ATG, 2 == TATGA, ...")
	output := flag.String("o", "stdout", "Output summary file.")
	btrim := flag.Int("bTrim", 0, "Trim ends of bed records in -b by #bp.")
	genomeCacheOutput := flag.String("genomeCacheOutput", "", "Output the results of genome context calculation to file to be used as input for future runs.")
	genomeCacheInput := flag.String("genomeCacheInput", "", "Input a genome cache file generated from a previous run to speed up execution.")
	verbose := flag.Int("v", 0, "Verbose output by setting to >0.")
	flag.Parse()

	if *input == "" || *ref == "" || *bedfile == "" {
		usage()
		log.Fatalln("ERROR: must have inputs for -i, -b, and -r")
	}

	burden.Burden(*input, *bedfile, *ref, *output, *pad, *btrim, *verbose, *genomeCacheInput, *genomeCacheOutput)
}
