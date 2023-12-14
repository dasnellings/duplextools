package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"path"
	"strings"
	"sync"
	"time"
)

func usage() {
	fmt.Print(
		"duplexMultiomeSplit - Split duplex multiome BAM file by cell barcode tag.\n" +
			"Usage:\n" +
			"annotateReadFamilies [options] -i input.bam -b barcodes.cells -strand1 s1.barcodes -strand2 s2.barcodes\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input bam file.")
	barcodesFile := flag.String("b", "", "Cell barcodes file.")
	strand1File := flag.String("strand1", "", "Strand 1 barcodes file. 1 barcode per line.")
	strand2File := flag.String("strand2", "", "Strand 2 barcodes file. 1 barcode per line.")
	outputDir := flag.String("outputDir", "barcode_split_bams", "Directory to output split bam files.")
	flag.Parse()

	if *input == "" || *barcodesFile == "" || *strand1File == "" || *strand2File == "" {
		usage()
		log.Fatal("ERROR: Must input a bam file and 3 barcodes files.")
	}

	duplexMultiomeSplit(*input, *barcodesFile, *outputDir, *strand1File, *strand2File)
}

func duplexMultiomeSplit(input, barcodesFile, outputDir, strand1File, strand2File string) {
	startTime := time.Now()
	wg := new(sync.WaitGroup)
	var err error
	err = os.Mkdir(outputDir, 0755)
	if err != nil {
		log.Fatalf("ERROR: output directory '%s' already exists.", outputDir)
	}
	records, header := sam.GoReadToChan(input)

	// make outputs writing to buffers in memory
	chanMap := make(map[string][2]chan<- sam.Sam)
	barcodes := readBarcodes(barcodesFile)
	var barcode string
	for _, barcode = range barcodes {
		wg.Add(2)
		chanMap[barcode] = [2]chan<- sam.Sam{
			newWriter(header, outputDir, path.Base(strings.TrimSuffix(input, ".bam")), barcode, wg, 1),
			newWriter(header, outputDir, path.Base(strings.TrimSuffix(input, ".bam")), barcode, wg, 2),
		}
	}

	strand1Barcodes := getStrandMap(strand1File)
	strand2Barcodes := getStrandMap(strand2File)

	var found bool
	var cellBarcode, strandBarcode string
	var value interface{}
	var currWriter chan<- sam.Sam
	var recordsProcessed int
	var unhandledRecords int
	for r := range records {
		recordsProcessed++
		if recordsProcessed%1000000 == 0 {
			log.Printf("Total Reads Processed: %d\n", recordsProcessed)
		}

		value, found, err = sam.QueryTag(r, "CB")
		exception.PanicOnErr(err)
		if !found {
			//log.Panicln("ERROR: cell barcode not found in following record\n", r)
			unhandledRecords++
			continue
		}
		cellBarcode = value.(string)

		value, found, err = sam.QueryTag(r, "BC")
		exception.PanicOnErr(err)
		if !found {
			//log.Panicln("ERROR: strand barcode not found in following record\n", r)
			unhandledRecords++
			continue
		}
		strandBarcode = value.(string)

		switch {
		case strand1Barcodes[strandBarcode]:
			currWriter = chanMap[cellBarcode][0]
		case strand2Barcodes[strandBarcode]:
			currWriter = chanMap[cellBarcode][1]
		}

		if currWriter == nil {
			unhandledRecords++
			continue
		}

		currWriter <- r
	}

	for _, c := range chanMap {
		close(c[0])
		close(c[1])
	}

	wg.Wait()
	elapsedTime := time.Since(startTime)
	log.Printf("\nProcessed Records: %d\nUnhandled Records: %d\nTime Elapsed: %.1fsec\n", recordsProcessed, unhandledRecords, elapsedTime.Seconds())
}

func readBarcodes(file string) []string {
	ans := make([]string, 0)
	f, err := os.Open(file)
	exception.PanicOnErr(err)
	defer f.Close()
	s := bufio.NewScanner(f)
	for s.Scan() {
		ans = append(ans, s.Text())
	}
	err = s.Err()
	exception.PanicOnErr(err)
	return ans
}

func writeBufferToFile(b *bytes.Buffer, dir, filename, barcode string, strand int) {
	var file *os.File
	var err error
	var filepath string
	if b.Len() == 0 {
		return
	}
	filepath = path.Join(dir, fmt.Sprintf("%s.%s.strand%d.bam", filename, barcode, strand))
	file, err = os.OpenFile(filepath, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0600)
	exception.PanicOnErr(err)
	_, err = file.Write(b.Bytes())
	exception.PanicOnErr(err)
	b.Reset()
	err = file.Close()
	exception.PanicOnErr(err)
}

func newWriter(header sam.Header, dir, filename, barcode string, wg *sync.WaitGroup, strand int) chan<- sam.Sam {
	in := make(chan sam.Sam, 100)
	go func(<-chan sam.Sam) {
		b := new(bytes.Buffer)
		w := sam.NewBamWriter(b, header)
		for r := range in {
			sam.WriteToBamFileHandle(w, r, 0)
			writeBufferToFile(b, dir, filename, barcode, strand)
		}
		err := w.Close()
		exception.PanicOnErr(err)
		writeBufferToFile(b, dir, filename, barcode, strand)
		wg.Done()
	}(in)
	return in
}

func getStrandMap(file string) map[string]bool {
	barcodes := readBarcodes(file)
	m := make(map[string]bool)
	for i := range barcodes {
		m[barcodes[i]] = true
	}
	return m
}
