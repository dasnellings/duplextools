package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/burden"
	"github.com/vertgenlab/gonomics/exception"
)

func burdenUsage(burdenFlags *flag.FlagSet) {
	fmt.Print(
		"burden - correct mutation burden for trinucleotide (or other) context bias and allele count\n" +
			"\tApplies method developed for NanoSeq (for more information see Abascal et al. 2021 PMID: 33911282)\n\n" +
			"Usage:\n" +
			"  duplextools burden [options] -i input.vcf -b readFamilies.bed -r reference.fasta > summary.tsv\n\n" +
			"Options:\n")
	burdenFlags.PrintDefaults()
}

func runBurden(args []string) {
	var err error
	burdenFlags := flag.NewFlagSet("burden", flag.ExitOnError)

	input := burdenFlags.String("i", "", "Input VCF file with final variant calls (somatic only, after filters).")
	ref := burdenFlags.String("r", "", "Reference FASTA file. Must be indexed (.fai).")
	bedfile := burdenFlags.String("b", "", "Bed file of read families used for calling. "+
		"In almost all cases, this should be the -b output of 'duplextools call'. If you are not using the file from 'call', "+
		"it is cricial that the bed regions used only include bases actually considered for variant calling. "+
		"If ends were trimmed from read families for variant calling, it should be reflected in the input bed regions or bTrim. "+
		"DO NOT INCLUDE FAMILIES THAT WERE EXCLUDED FROM ANALYSIS!!!")
	pad := burdenFlags.Int("pad", 1, "Number of context bases on either side of variant (e.g. 0 == T, 1 == ATA, 2 == AATAA, ...)")
	output := burdenFlags.String("o", "stdout", "Output summary file.")
	btrim := burdenFlags.Int("bTrim", 0, "Trim ends of bed records in -b by #bp.")
	genomeCacheOutput := burdenFlags.String("genomeCacheOutput", "", "Output the results of genome context calculation to file to be used as input for future runs.")
	genomeCacheInput := burdenFlags.String("genomeCacheInput", "", "Input a genome cache file generated from a previous run to speed up execution.")
	verbose := burdenFlags.Int("v", 0, "Verbose output by setting to >0.")

	err = burdenFlags.Parse(args)
	exception.PanicOnErr(err)
	burdenFlags.Usage = func() { burdenUsage(burdenFlags) }

	if *input == "" || *ref == "" || *bedfile == "" {
		burdenFlags.Usage()
		errExit("\nERROR: must have inputs for -i, -b, and -r")
	}

	burden.Burden(*input, *bedfile, *ref, *output, *pad, *btrim, *verbose, *genomeCacheInput, *genomeCacheOutput)
}
