package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/pair"
	"github.com/vertgenlab/gonomics/exception"
)

func pairUsage(pairFlags *flag.FlagSet) {
	fmt.Print(
		"pair - identify read families by duplex barcodes and start/end site and record family in RF tag for each read\n\n" +
			"Usage:\n" +
			"  duplextools pair [options] -i input.bam > output.bam\n\n" +
			"Options:\n")
	pairFlags.PrintDefaults()
}

func runPair(args []string) {
	var err error
	pairFlags := flag.NewFlagSet("pair", flag.ExitOnError)

	input := pairFlags.String("i", "", "Input bam file. Must be coordinate sorted.")
	output := pairFlags.String("o", "stdout", "Output bam file.")
	bed := pairFlags.String("bed", "", "Output a bed file with the region covered by each read family. May significantly increase memory usage.")
	strict := pairFlags.Bool("strict", false, "Require perfect barcode match for family inclusion. Disables position matching. Use for high density data.")
	tolerance := pairFlags.Int("tolerance", 50, "Deviation from exact start match to be considered for inclusion in read family. 0 means perfect match. Low values are best for dense data, and high values are best for sparse data.")
	strictPosMatching := pairFlags.Bool("strictPosMatching", false, "For a read to be included in a read family, the start of both reads in a pair must exactly match the read family.")
	minMapQ := pairFlags.Int("minMapQ", 20, "Minimum mapping quality.")

	err = pairFlags.Parse(args)
	exception.PanicOnErr(err)
	pairFlags.Usage = func() { pairUsage(pairFlags) }

	if *input == "" {
		pairFlags.Usage()
		errExit("\nERROR: must input a coordinate sorted bam file with -i")
	}

	pair.Pair(*input, *output, *tolerance, *strict, *strictPosMatching, *bed, uint8(*minMapQ))
}
