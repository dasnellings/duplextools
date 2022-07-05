package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/families"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
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
	strict := flag.Bool("strict", false, "Require perfect barcode match for family inclusion. Disables position matching. Use for high density data.")
	tolerance := flag.Int("t", 1000, "Deviation from exact start match to be considered for inclusion in read family. 0 means perfect match. Low values are best for dense data, and high values are best for sparse data.")
	flag.Parse()

	if *input == "" {
		usage()
		log.Fatal("ERROR: Must input a coordinate sorted bam file.")
	}

	annotateReadFamilies(*input, *output, *tolerance, *strict)
}

func annotateReadFamilies(input, output string, tolerance int, strict bool) {
	reads, header := sam.GoReadToChan(input)
	if header.Metadata.SortOrder[0] != sam.Coordinate {
		log.Fatal("ERROR: Input file must be coordinate sorted.")
	}

	reads = families.GoAnnotate(reads, tolerance, !strict)

	out := fileio.EasyCreate(output)
	bw := sam.NewBamWriter(out, header)

	for r := range reads {
		if r.RName == "" {
			continue
		}
		sam.WriteToBamFileHandle(bw, r, 0)
	}

	err := bw.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}
