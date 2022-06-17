package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

// SAM format uses ascii offset of 33 to make everything start with individual characters
// without adding 33 you get values like spaces and newlines
const asciiOffset uint8 = 33

func usage() {
	fmt.Print(
		"mcsFqToBam - Process raw FASTQ files generated with META-CS to an unmapped BAM with barcode tags.\n\n" +
			"Usage:\n" +
			"  mcsFqToBam [options] -1 r1.fq.gz -2 r2.fq.gz -o output.bam\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func fqToSam(fq *fastq.Fastq, sam *sam.Sam, firstInPair bool) {
	sam.QName = fq.Name
	sam.Seq = fq.Seq
	for i := range fq.Qual {
		fq.Qual[i] += asciiOffset
	}
	sam.Qual = string(fq.Qual)
	if firstInPair {
		sam.Flag = 77
	} else {
		sam.Flag = 141
	}
}

func mcsFqToBam(r1File, r2File, outFile string) {
	readPairs := make(chan fastq.PairedEnd, 1000)
	go fastq.PairedEndToChan(r1File, r2File, readPairs)

	o := fileio.EasyCreate(outFile)
	bw := sam.NewBamWriter(o, sam.GenerateHeader(nil, nil, sam.Unsorted, sam.None))

	var pair fastq.PairedEnd
	var s1, s2 sam.Sam
	for pair = range readPairs {
		fqToSam(&pair.Fwd, &s1, true)
		fqToSam(&pair.Rev, &s2, false)

		sam.WriteToBamFileHandle(bw, s1, 0)
		sam.WriteToBamFileHandle(bw, s2, 0)
	}

	err := bw.Close()
	exception.PanicOnErr(err)
	err = o.Close()
	exception.PanicOnErr(err)
}

func main() {
	r1 := flag.String("1", "", "FASTQ file containing R1 reads. May be gzipped.")
	r2 := flag.String("2", "", "FASTQ file containing R2 reads. May be gzipped.")
	outfile := flag.String("o", "stdout", "Output BAM file.")
	flag.Parse()
	flag.Usage = usage

	if *r1 == "" || *r2 == "" {
		flag.Usage()
		log.Fatalln("ERROR: must input read1 file and read2 file")
	}

	mcsFqToBam(*r1, *r2, *outfile)
}
