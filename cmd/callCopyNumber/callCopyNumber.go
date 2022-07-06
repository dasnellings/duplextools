package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"math/rand"
	"sort"
	"strings"
)

const debug int = 1

var countMat []int

func usage() {
	fmt.Print(
		"callCopyNumber - Call copy number from read family collapsed BED file generated with annotateReadFamilies.\n" +
			"Usage:\n" +
			"callCopyNumber [options] -i readFamilies.bed\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input bed file. Must be coordinate sorted.")
	output := flag.String("o", "stdout", "Output bedgraph file.")
	minOverlap := flag.Int("minOverlap", 10, "Minimum overlap of two read families to be considered from different alleles. Tn5 generates a 9bp overlap.")
	maxBedLength := flag.Int("maxBedLen", 2000, "Maximum size of a bed record for inclusion in analysis.")
	minReads := flag.Int("minReads", 3, "Minimum size of read family for inclusion in analysis.")
	mergeIdenticalPos := flag.Bool("merge", true, "Merge bed records with identical starts OR identical ends.")
	//minReadsPerFamily := flag.Int("minReads", 1, "Minimum number of reads in a read family for inclusion in analysis.")
	flag.Parse()

	if *input == "" {
		usage()
		log.Fatal("ERROR: Must input a coordinate sorted bed file.")
	}

	callCopyNumber(*input, *output, *minOverlap-1, *maxBedLength, *minReads, *mergeIdenticalPos)
}

func callCopyNumber(input, output string, trimLen, maxBedLen, minReads int, merge bool) {
	var err error
	var prevChrom string
	var prevStart int
	var overlapSet []bed.Bed
	records := bed.GoReadToChan(input)
	out := fileio.EasyCreate(output)

	var debugOut io.WriteCloser
	if debug == 1 {
		debugOut = fileio.EasyCreate("debug_" + strings.TrimRight(output, ".bedgraph") + ".bed")
	}
	for b := range records {
		if b.Chrom == prevChrom && b.ChromStart < prevStart {
			log.Fatalln("ERROR: Input bed file is not coordinate sorted.", b, prevStart)
		}
		prevChrom = b.Chrom
		prevStart = b.ChromStart

		err = trim(&b, trimLen)
		if err != nil || b.ChromEnd-b.ChromStart > maxBedLen || b.Score < minReads {
			continue
		}

		switch {
		case len(overlapSet) == 0: // nothing in overlap set
			overlapSet = append(overlapSet, b)

		case !anyOverlaps(overlapSet, b): // no overlaps with set
			writeOverlapsCounting(out, overlapSet, merge, debugOut)
			overlapSet = overlapSet[:0]
			overlapSet = append(overlapSet, b)

		default: // overlaps the set
			overlapSet = append(overlapSet, b)
		}
	}
	writeOverlapsCounting(out, overlapSet, merge, debugOut)
	err = out.Close()
	exception.PanicOnErr(err)
	if debug == 1 {
		err = debugOut.Close()
		exception.PanicOnErr(err)
	}
}

func writeOverlaps(out io.Writer, set []bed.Bed, debugOut io.Writer) {
	t := rand.Uint32()
	if len(set) == 0 {
		return
	}

	if debug == 1 {
		for i := range set {
			fmt.Fprintf(debugOut, "%s\t%d\t%d\t%d\n", set[i].Chrom, set[i].ChromStart, set[i].ChromEnd, t)
		}
	}
	var curr bed.Bed
	var currCount int = 1
	curr.Chrom = set[0].Chrom
	curr.ChromStart = set[0].ChromStart
	curr.ChromEnd = set[0].ChromEnd
	for i := 1; i < len(set); i++ {
		if set[i].ChromStart > curr.ChromStart {
			if set[i].ChromStart < curr.ChromEnd {
				curr.ChromEnd = set[i].ChromStart
			}
			fmt.Fprintf(out, "%s\t%d\t%d\t%d\n", curr.Chrom, curr.ChromStart, curr.ChromEnd, currCount)
			i = 1
			currCount = 1
			resetStack(&set, &curr)
		} else {
			if set[i].ChromEnd < curr.ChromEnd {
				curr.ChromEnd = set[i].ChromEnd
			}
			currCount++
		}
	}
}

func writeOverlapsCounting(out io.Writer, set []bed.Bed, merge bool, debugOut io.Writer) {
	t := rand.Uint32()
	if len(set) == 0 {
		return
	}

	if merge {
		set = mergeIdenticalEndpoints(set)
	}

	if debug == 1 {
		for i := range set {
			fmt.Fprintf(debugOut, "%s\t%d\t%d\t%d\n", set[i].Chrom, set[i].ChromStart, set[i].ChromEnd, t)
		}
	}

	var start, end, currStart, currEnd int
	start = set[0].ChromStart
	for i := range set {
		if set[i].ChromEnd > end {
			end = set[i].ChromEnd
		}
	}
	length := end - start

	if cap(countMat) < length {
		countMat = make([]int, length)
	} else {
		countMat = countMat[0:length]
		for i := range countMat {
			countMat[i] = 0
		}
	}

	var i, j int
	for i = range set {
		currStart = set[i].ChromStart - start
		currEnd = set[i].ChromEnd - start

		for j = currStart; j < currEnd; j++ {
			countMat[j]++
		}
	}

	var curr bed.Bed
	curr.Chrom = set[0].Chrom
	curr.ChromStart = set[0].ChromStart
	for i = 1; i < len(countMat); i++ {
		if countMat[i] != countMat[i-1] {
			curr.ChromEnd = start + i
			fmt.Fprintf(out, "%s\t%d\t%d\t%d\n", curr.Chrom, curr.ChromStart, curr.ChromEnd, countMat[i-1])
			curr.ChromStart = curr.ChromEnd
		}
	}
	curr.ChromEnd = end
	fmt.Fprintf(out, "%s\t%d\t%d\t%d\n", curr.Chrom, curr.ChromStart, curr.ChromEnd, countMat[len(countMat)-1])
}

func anyOverlaps(set []bed.Bed, query bed.Bed) bool {
	for i := range set {
		if bed.Overlap(set[i], query) {
			return true
		}
	}
	return false
}

func trim(b *bed.Bed, len int) error {
	b.ChromStart += len
	b.ChromEnd -= len
	if b.ChromStart > b.ChromEnd {
		return errors.New("start gth end")
	}
	return nil
}

func resetStack(set *[]bed.Bed, curr *bed.Bed) {
	for i := 0; i < len(*set); i++ {
		switch {
		case (*set)[i].ChromEnd <= curr.ChromEnd:
			*set = slices.Delete(*set, i, i+1)
		case (*set)[i].ChromStart < curr.ChromEnd:
			(*set)[i].ChromStart = curr.ChromEnd
		}
	}
	curr.ChromStart = (*set)[0].ChromStart
	curr.ChromEnd = (*set)[0].ChromEnd
}

func mergeIdenticalEndpoints(set []bed.Bed) []bed.Bed {
	for i := 1; i < len(set); i++ {
		if set[i].ChromStart == set[i-1].ChromStart {
			set[i].ChromEnd = max(set[i].ChromEnd, set[i-1].ChromEnd)
			set = slices.Delete(set, i-1, i)
		}
	}

	sort.Slice(set, func(i, j int) bool {
		return set[i].ChromEnd < set[j].ChromEnd
	})

	for i := 1; i < len(set); i++ {
		if set[i].ChromEnd == set[i-1].ChromEnd {
			set[i].ChromStart = min(set[i].ChromStart, set[i-1].ChromStart)
			set = slices.Delete(set, i-1, i)
		}
	}

	sort.Slice(set, func(i, j int) bool {
		return set[i].ChromStart < set[j].ChromStart
	})

	return set
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
