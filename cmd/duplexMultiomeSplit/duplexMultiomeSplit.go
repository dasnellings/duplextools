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
	"time"
)

func usage() {
	fmt.Print(
		"duplexMultiomeSplit - Split duplex multiome BAM file by cell barcode tag.\n" +
			"Usage:\n" +
			"annotateReadFamilies [options] -i input.bam -b barcodes.cells\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input bam file.")
	barcodesFile := flag.String("b", "", "Barcodes file.")
	outputDir := flag.String("outputDir", "barcode_split_bams", "Directory to output split bam files.")
	flag.Parse()

	if *input == "" || *barcodesFile == "" {
		usage()
		log.Fatal("ERROR: Must input a bam file and a barcodes file.")
	}

	duplexMultiomeSplit(*input, *barcodesFile, *outputDir)
}

func duplexMultiomeSplit(input, barcodesFile, outputDir string) {
	startTime := time.Now()
	var err error
	os.Mkdir(outputDir, 0755)
	records, header := sam.GoReadToChan(input)

	// make outputs writing to buffers in memory
	barcodeMap := make(map[string][2]*sam.BamWriter)
	writerMap := make(map[string][2]*bytes.Buffer)
	barcodes := readBarcodes(barcodesFile)
	var barcode string
	for _, barcode = range barcodes {
		writerMap[barcode] = [2]*bytes.Buffer{new(bytes.Buffer), new(bytes.Buffer)}
		barcodeMap[barcode] = [2]*sam.BamWriter{sam.NewBamWriter(writerMap[barcode][0], header), sam.NewBamWriter(writerMap[barcode][1], header)}
	}

	var found bool
	var cellBarcode, strandBarcode string
	var value interface{}
	var currWriter *sam.BamWriter
	var recordsProcessed int
	var unhandledRecords int
	for r := range records {
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

		switch strandBarcode {
		case "AAACGGCG", "CCTACCAT", "GGCGTTTC", "TTGTAAGA":
			currWriter = barcodeMap[cellBarcode][0]
		case "AGGCTACC", "CTAGCTGT", "GCCAACAA", "TATTGGTG":
			currWriter = barcodeMap[cellBarcode][1]
		}

		if currWriter == nil {
			unhandledRecords++
			continue
		}

		sam.WriteToBamFileHandle(currWriter, r, 0)
		recordsProcessed++

		if recordsProcessed%100000 == 0 {
			log.Printf("Total Reads Processed: %d\n", recordsProcessed)
			appendBuffersToFile(writerMap, outputDir, strings.TrimSuffix(input, ".bam"))
		}
	}

	for _, writer := range barcodeMap {
		err = writer[0].Close()
		exception.PanicOnErr(err)
		err = writer[1].Close()
		exception.PanicOnErr(err)
	}
	appendBuffersToFile(writerMap, outputDir, strings.TrimSuffix(input, ".bam"))

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

func appendBuffersToFile(writerMap map[string][2]*bytes.Buffer, dir, filename string) {
	var file *os.File
	var err error
	var filepath string
	var i int
	for barcode, writers := range writerMap {
		for i = range writers {
			if writers[i].Len() == 0 {
				continue
			}
			filepath = path.Join(dir, fmt.Sprintf("%s.%s.strand%d.bam", filename, barcode, i+1))
			file, err = os.OpenFile(filepath, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0600)
			exception.PanicOnErr(err)
			_, err = file.Write(writers[i].Bytes())
			exception.PanicOnErr(err)
			writers[i].Reset()
			err = file.Close()
			exception.PanicOnErr(err)
		}
	}
}
