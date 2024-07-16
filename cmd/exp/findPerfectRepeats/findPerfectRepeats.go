package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/repeats"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"findPerfectRepeats - Parse TRF output to identify perfect repeats.\n" +
			"Usage:\n" +
			"rmTag [options] -i trf_output.txt -r reference.fasta > output.bed\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input Tandem Repeat Finder text file.")
	output := flag.String("o", "stdout", "Output `BED` file.")
	ref := flag.String("r", "", "Reference `FASTA` file, must be indexed.")
	minRepeatUnits := flag.Int("minRepeatUnits", 10, "Minimum number of repeated units to be included in output.")
	minUnitLen := flag.Int("minUnitLen", 1, "Minimum length of repeat unit to be included in output.")
	maxUnitLen := flag.Int("maxUnitLen", 10, "Maximum length of repeat unit to be included in output.")
	maxTotalLen := flag.Int("maxTotalLen", 75, "Maximum total length of repeat.")
	distToUnmasked := flag.Int("maxDistToUnmasked", 20, "Maximum distance from both ends of repeat to unmasked sequence (as determined by case in -r) to be included in output. -1 to disable")
	flag.Parse()

	if *input == "" || *ref == "" {
		usage()
		log.Fatal("ERROR: Must have values for -i and -r.")
	}

	findPerfectRepeats(*input, *output, *ref, *minRepeatUnits, *minUnitLen, *maxUnitLen, *maxTotalLen, *distToUnmasked)
}

func findPerfectRepeats(input, output, reference string, minRepeatUnits, minUnitLen, maxUnitLen, maxTotalLen, distToUnmasked int) {
	records := repeats.GoReadToChan(input)
	out := fileio.EasyCreate(output)
	defer cleanup(out)
	ref := fasta.NewSeeker(reference, "")
	defer cleanup(ref)

	var curr bed.Bed
	var numRepeats int
	var repeatSeq, refSeq []dna.Base
	for r := range records {
		if strings.ContainsAny(r.Chr, "_*v") {
			continue
		}
		curr, numRepeats, repeatSeq = repeats.FindPerfectRepeat(ref, &r)
		if numRepeats < minRepeatUnits || len(repeatSeq) < minUnitLen || len(repeatSeq) > maxUnitLen || len(repeatSeq)*numRepeats > maxTotalLen {
			continue
		}
		//if float64(r.End-r.Start) > 1.25*float64(len(repeatSeq)*numRepeats) {
		//	continue
		//}

		if distToUnmasked > -1 {
			refSeq, _ = fasta.SeekByName(ref, r.Chr, curr.ChromStart-distToUnmasked, curr.ChromStart)
			if !containsUnmasked(refSeq) {
				continue
			}
			refSeq, _ = fasta.SeekByName(ref, r.Chr, curr.ChromEnd, curr.ChromEnd+distToUnmasked)
			if !containsUnmasked(refSeq) {
				continue
			}
		}
		bed.WriteBed(out, curr)
	}
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}

func containsUnmasked(seq []dna.Base) bool {
	for i := range seq {
		switch seq[i] {
		case dna.A, dna.C, dna.G, dna.T:
			return true
		}
	}
	return false
}
