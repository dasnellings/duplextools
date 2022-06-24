package main

import (
	"fmt"
	"github.com/dasnellings/MCS_MS/barcode"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
)

func main() {
	infile := os.Args[1]
	reads, header := sam.GoReadToChan(infile)
	if header.Metadata.SortOrder[0] != sam.Coordinate {
		log.Fatal("ERROR: input file must be coordinate sorted")
	}

	var currChrom string
	var currStart uint32
	var totalReads, duplexSites, totalSites int
	var bcFor, bcRev, currBcFor, currBcRev string
	for r := range reads {
		totalReads++
		if currStart == r.Pos && currChrom == r.RName {
			currBcFor, currBcRev = barcode.Get(r)
			if currBcFor == bcRev && currBcRev == bcFor {
				duplexSites++
				bcFor = "dummy"
				bcRev = "dummy"
			}
		} else {
			totalSites++
			currStart = r.Pos
			currChrom = r.RName
			bcFor, bcRev = barcode.Get(r)
		}
	}

	fmt.Printf("Total Reads:\t%d\n", totalReads)
	fmt.Printf("Total Sites:\t%d\n", totalSites)
	fmt.Printf("Duplex Sites:\t%d\n", duplexSites)
	fmt.Printf("Duplex Fraction:\t%f\n", float64(duplexSites)/float64(totalSites))
}
