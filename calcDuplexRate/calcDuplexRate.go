package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/barcode"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

type output struct {
	chrom    string
	totReads int
	totSites int
	dupSites int
}

const startTolerance int = 0

func main() {
	tolerance := flag.Int("t", 0, "Deviation from exact start match to be considered same allele. 0 means perfect match.")
	infile := flag.String("i", "", "Input coordinate sorted BAM or SAM file.")
	update := flag.Int("u", 0, "Print duplex rate in chunks, ever INT reads. 0 only reports after all data is read.")
	flag.Parse()

	if *infile == "" {
		flag.PrintDefaults()
		log.Fatal("ERROR: must input coordinate sorted BAM or SAM file")
	}
	reads, header := sam.GoReadToChan(*infile)
	if header.Metadata.SortOrder[0] != sam.Coordinate {
		log.Fatal("ERROR: input file must be coordinate sorted")
	}
	updateFreq := *update
	var chunkStartChrom, chunkEndChrom string
	var chunkTotalSites, chunkDuplexSites int
	var chunkStart, chunkEnd uint32
	if updateFreq > 0 {
		fmt.Println("Region\tTotalSites\tDuplexSites\tDuplexFraction")
	}

	var a []output

	startTolerance := uint32(*tolerance)
	var currChrom string
	var currStart uint32
	var totalReads, duplexSites, totalSites int
	var bcFor, bcRev, currBcFor, currBcRev string
	for r := range reads {

		if r.RName != currChrom {
			a = append(a, output{currChrom, totalReads, totalSites, duplexSites})
			currChrom = r.RName
			totalReads = 0
			totalSites = 0
			duplexSites = 0
		}

		totalReads++
		if (currStart >= r.Pos-startTolerance && currStart <= r.Pos+startTolerance) && currChrom == r.RName {
			currBcFor, currBcRev = barcode.Get(r)
			if currBcFor == bcRev && currBcRev == bcFor {
				duplexSites++
				chunkDuplexSites++
				bcFor = "dummy"
				bcRev = "dummy"
			}
		} else {
			totalSites++
			chunkTotalSites++
			currStart = r.Pos
			currChrom = r.RName
			bcFor, bcRev = barcode.Get(r)
		}

		if updateFreq > 0 && totalReads%updateFreq == 0 {
			chunkEndChrom = r.RName
			chunkEnd = r.Pos
			if chunkStartChrom == "" {
				chunkStartChrom = chunkEndChrom
			}
			fmt.Printf("%s:%d-%s:%d\t%d\t%d\t%f\n", chunkStartChrom, chunkStart, chunkEndChrom, chunkEnd, chunkTotalSites, chunkDuplexSites, float64(chunkDuplexSites)/float64(chunkTotalSites))
			chunkStartChrom = chunkEndChrom
			chunkStart = chunkEnd
			chunkTotalSites = 0
			chunkDuplexSites = 1
		}
	}

	fmt.Println("Chromosome\tTotal Reads\tTotal Sites\tDuplex Fraction")
	for i := range a {
		if a[i].chrom == "" {
			continue
		}
		fmt.Printf("%s\t%d\t%d\t%d\t%f\n", a[i].chrom, a[i].totReads, a[i].totSites, a[i].totReads, float64(a[i].dupSites)/float64(a[i].totSites))
	}
}
