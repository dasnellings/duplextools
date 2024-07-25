package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/pair"
	"log"
)

func usage() {
	fmt.Print(
		"annotateReadFamilies - Collapse read families from META-CS data and record families in RF tag for each read.\n" +
			"Usage:\n" +
			"annotateReadFamilies [options] -i input.bam > output.bam\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input bam file. Must be coordinate sorted.")
	output := flag.String("o", "stdout", "Output bam file.")
	bed := flag.String("bed", "", "Output a bed file with the region covered by each read family. May significantly increase memory usage.")
	strict := flag.Bool("strict", false, "Require perfect barcode match for family inclusion. Disables position matching. Use for high density data.")
	tolerance := flag.Int("tolerance", 50, "Deviation from exact start match to be considered for inclusion in read family. 0 means perfect match. Low values are best for dense data, and high values are best for sparse data.")
	strictPosMatching := flag.Bool("strictPosMatching", false, "For a read to be included in a read family, the start of both reads in a pair must exactly match the read family.")
	minMapQ := flag.Int("minMapQ", 20, "Minimum mapping quality.")
	flag.Parse()

	if *input == "" {
		usage()
		log.Fatal("ERROR: Must input a coordinate sorted bam file.")
	}

	pair.Pair(*input, *output, *tolerance, *strict, *strictPosMatching, *bed, uint8(*minMapQ))
}
