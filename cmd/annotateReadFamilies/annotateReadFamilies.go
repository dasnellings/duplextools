package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/families"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"log"
	"strings"
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
	tolerance := flag.Int("t", 1000, "Deviation from exact start match to be considered for inclusion in read family. 0 means perfect match. Low values are best for dense data, and high values are best for sparse data.")
	flag.Parse()

	if *input == "" {
		usage()
		log.Fatal("ERROR: Must input a coordinate sorted bam file.")
	}

	annotateReadFamilies(*input, *output, *tolerance, *strict, *bed)
}

type minimalBed struct {
	chr    string
	start  int
	end    int
	family string
}

func annotateReadFamilies(input, output string, tolerance int, strict bool, bed string) {
	var err error
	reads, header := sam.GoReadToChan(input)
	if header.Metadata.SortOrder[0] != sam.Coordinate {
		log.Fatal("ERROR: Input file must be coordinate sorted.")
	}

	reads = families.GoAnnotate(reads, tolerance, !strict)

	out := fileio.EasyCreate(output)
	bw := sam.NewBamWriter(out, header)

	var bedOut io.WriteCloser
	m := make(map[string]minimalBed)
	var rf string
	var mb minimalBed
	if bed != "" {
		bedOut = fileio.EasyCreate(bed)
	}

	for r := range reads {
		if r.RName == "" {
			continue
		}
		sam.WriteToBamFileHandle(bw, r, 0)

		if bed != "" {
			rf = getRF(&r)
			mb = m[rf]
			mb.chr = r.RName
			mb.family = rf
			if mb.start == 0 || mb.start > r.GetChromStart() {
				mb.start = r.GetChromStart()
			}
			if mb.end < r.GetChromEnd() {
				mb.end = r.GetChromEnd()
			}
			m[rf] = mb
		}
	}

	if bed != "" {
		for _, b := range m {
			fmt.Fprintf(bedOut, "%s\t%d\t%d\t%s", b.chr, b.start, b.end, b.family)
		}
		err = bedOut.Close()
		exception.PanicOnErr(err)
	}

	err = bw.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func getRF(r *sam.Sam) string {
	idx := strings.Index(r.Extra, "RF:Z:")
	if idx == -1 {
		return ""
	}
	return r.Extra[idx+5:]
}
