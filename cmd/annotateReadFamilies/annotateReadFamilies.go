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
	"sort"
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
	tolerance := flag.Int("tolerance", 50, "Deviation from exact start match to be considered for inclusion in read family. 0 means perfect match. Low values are best for dense data, and high values are best for sparse data.")
	minMapQ := flag.Int("minMapQ", 20, "Minimum mapping quality.")
	flag.Parse()

	if *input == "" {
		usage()
		log.Fatal("ERROR: Must input a coordinate sorted bam file.")
	}

	annotateReadFamilies(*input, *output, *tolerance, *strict, *bed, uint8(*minMapQ))
}

type minimalBed struct {
	chr    string
	start  int
	end    int
	family string
	count  int
}

func annotateReadFamilies(input, output string, tolerance int, strict bool, bed string, minMapQ uint8) {
	var err error
	reads, header := sam.GoReadToChan(input)
	if header.Metadata.SortOrder[0] != sam.Coordinate {
		log.Fatal("ERROR: Input file must be coordinate sorted.")
	}

	reads = families.GoAnnotate(reads, tolerance, !strict)

	out := fileio.EasyCreate(output)
	bw := sam.NewBamWriter(out, header)

	var bedOut io.WriteCloser
	m := make(map[string]*minimalBed)
	var rf string
	var mb *minimalBed
	if bed != "" {
		bedOut = fileio.EasyCreate(bed)
	}
	var prevChrom string
	var readCount int
	var bedToWrite []*minimalBed

	for r := range reads {
		if r.RName == "" || r.MapQ < minMapQ {
			continue
		}
		sam.WriteToBamFileHandle(bw, r, 0)

		if bed == "" {
			continue
		}
		readCount++
		if r.RName != prevChrom {
			for k, b := range m {
				bedToWrite = append(bedToWrite, b)
				delete(m, k)
			}
			sort.Slice(bedToWrite, func(i, j int) bool {
				switch {
				case bedToWrite[i].chr < bedToWrite[j].chr:
					return true
				case bedToWrite[i].chr > bedToWrite[j].chr:
					return false
				case bedToWrite[i].start < bedToWrite[j].start:
					return true
				case bedToWrite[i].start > bedToWrite[j].start:
					return false
				case bedToWrite[i].end < bedToWrite[j].end:
					return true
				default:
					return false
				}
			})
			for _, b := range bedToWrite {
				fmt.Fprintf(bedOut, "%s\t%d\t%d\t%s\t%d\n", b.chr, b.start, b.end, b.family, b.count)
			}
			bedToWrite = bedToWrite[:0]
		}

		rf = getRF(&r)
		mb = m[rf]
		if mb == nil {
			mb = new(minimalBed)
			m[rf] = mb
			mb.chr = r.RName
			mb.family = rf
		}
		if mb.start == 0 || mb.start > r.GetChromStart() {
			mb.start = r.GetChromStart()
		}
		if mb.end < r.GetChromEnd() {
			mb.end = r.GetChromEnd()
		}
		mb.count++
		prevChrom = r.RName

		if readCount%10000 == 0 { // write every 10000 reads
			for k, b := range m {
				if b.end < r.GetChromStart()-10000 { // only write if family is at least 10kb away
					bedToWrite = append(bedToWrite, b)
					delete(m, k)
				}
			}
			sort.Slice(bedToWrite, func(i, j int) bool {
				switch {
				case bedToWrite[i].chr < bedToWrite[j].chr:
					return true
				case bedToWrite[i].chr > bedToWrite[j].chr:
					return false
				case bedToWrite[i].start < bedToWrite[j].start:
					return true
				case bedToWrite[i].start > bedToWrite[j].start:
					return false
				case bedToWrite[i].end < bedToWrite[j].end:
					return true
				default:
					return false
				}
			})
			for _, b := range bedToWrite {
				fmt.Fprintf(bedOut, "%s\t%d\t%d\t%s\t%d\n", b.chr, b.start, b.end, b.family, b.count)
			}
			bedToWrite = bedToWrite[:0]
		}
	}

	if bed != "" {
		for k, b := range m {
			bedToWrite = append(bedToWrite, b)
			delete(m, k)
		}
		sort.Slice(bedToWrite, func(i, j int) bool {
			switch {
			case bedToWrite[i].chr < bedToWrite[j].chr:
				return true
			case bedToWrite[i].chr > bedToWrite[j].chr:
				return false
			case bedToWrite[i].start < bedToWrite[j].start:
				return true
			case bedToWrite[i].start > bedToWrite[j].start:
				return false
			case bedToWrite[i].end < bedToWrite[j].end:
				return true
			default:
				return false
			}
		})
		for _, b := range bedToWrite {
			fmt.Fprintf(bedOut, "%s\t%d\t%d\t%s\t%d\n", b.chr, b.start, b.end, b.family, b.count)
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
