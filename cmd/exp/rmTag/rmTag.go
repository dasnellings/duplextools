package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"rmTag - Remove a tag from read in a bam file.\n" +
			"Usage:\n" +
			"rmTag [options] -i input.bam -tag RF > output.bam\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input bam file.")
	output := flag.String("o", "stdout", "Output bam file.")
	tag := flag.String("tag", "", "Tag to remove")
	flag.Parse()

	if *input == "" || *tag == "" {
		usage()
		log.Fatal("ERROR: Must input a coordinate sorted bam file.")
	}

	t := *tag

	reads, header := sam.GoReadToChan(*input)
	out := fileio.EasyCreate(*output)
	bw := sam.NewBamWriter(out, header)
	var start, end, i int
	for r := range reads {
		sam.ParseExtra(&r)
		start = strings.Index(r.Extra, t+":")
		end = -1
		for i = start; i < len(r.Extra); i++ {
			if r.Extra[i] == '\t' {
				end = i
				break
			}
		}

		if end == -1 {
			r.Extra = r.Extra[:start-1]
		} else {
			r.Extra = r.Extra[:start-1] + r.Extra[end:]
		}

		sam.WriteToBamFileHandle(bw, r, 0)
	}

	err := bw.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}
