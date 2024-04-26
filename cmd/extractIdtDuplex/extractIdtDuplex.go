package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"os"
	"strings"
)

const (
	i7FivePrimeFlankFwd  string = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	i5FivePrimeFlankFwd  string = "AATGATACGGCGACCACCGAGATCTACAC"
	i7ThreePrimeFlankFwd string = "ATCTCGTATGCCGTCTTCTGCTTG"
	i5ThreePrimeFlankFwd string = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"

	i7ThreePrimeFlankRev string = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
	i5ThreePrimeFlankRev string = "GTGTAGATCTCGGTGGTCGCCGTATCATT"
	i7FivePrimeFlankRev  string = "CAAGCAGAAGACGGCATACGAGAT"
	i5FivePrimeFlankRev  string = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

	fwdStart string = "AATGATACGGCGACCACCGAGATCTACAC"
	revStart string = "CAAGCAGAAGACGGCATACGAGAT"
)

const (
	UMI_1  string = "GAGACGAT"
	UMI_2  string = "TTCCAAGG"
	UMI_3  string = "CGCATGAT"
	UMI_4  string = "ACGGAACA"
	UMI_5  string = "CGGCTAAT"
	UMI_6  string = "GCTATCCT"
	UMI_7  string = "TGGACTCT"
	UMI_8  string = "ATCCAGAG"
	UMI_9  string = "CTTAGGAC"
	UMI_10 string = "GTGCCATA"
	UMI_11 string = "TCGCTGTT"
	UMI_12 string = "TTCGTTGG"
	UMI_13 string = "AAGCACTG"
	UMI_14 string = "GTCGAAGA"
	UMI_15 string = "ACCACGAT"
	UMI_16 string = "GATTACCG"
	UMI_17 string = "GCACAACT"
	UMI_18 string = "GCGTCATT"
	UMI_19 string = "GAAGGAAG"
	UMI_20 string = "ACTGAGGT"
	UMI_21 string = "TGAAGACG"
	UMI_22 string = "GTTACGCA"
	UMI_23 string = "AGCGTGTT"
	UMI_24 string = "GATCGAGT"
	UMI_25 string = "TTGCGAAG"
	UMI_26 string = "CTGTTGAC"
	UMI_27 string = "GATGTGTG"
	UMI_28 string = "ACGTTCAG"
	UMI_29 string = "TTGCAGAC"
	UMI_30 string = "CAATGTGG"
	UMI_32 string = "ACTAGGAG"
)

var umiMap = map[string]string{
	"GAGACGAT": "UMI_1", "TTCCAAGG": "UMI_2", "CGCATGAT": "UMI_3",
	"ACGGAACA": "UMI_4", "CGGCTAAT": "UMI_5", "GCTATCCT": "UMI_6",
	"TGGACTCT": "UMI_7", "ATCCAGAG": "UMI_8", "CTTAGGAC": "UMI_9",
	"GTGCCATA": "UMI_10", "TCGCTGTT": "UMI_11", "TTCGTTGG": "UMI_12",
	"AAGCACTG": "UMI_13", "GTCGAAGA": "UMI_14", "ACCACGAT": "UMI_15",
	"GATTACCG": "UMI_16", "GCACAACT": "UMI_17", "GCGTCATT": "UMI_18",
	"GAAGGAAG": "UMI_19", "ACTGAGGT": "UMI_20", "TGAAGACG": "UMI_21",
	"GTTACGCA": "UMI_22", "AGCGTGTT": "UMI_23", "GATCGAGT": "UMI_24",
	"TTGCGAAG": "UMI_25", "CTGTTGAC": "UMI_26", "GATGTGTG": "UMI_27",
	"ACGTTCAG": "UMI_28", "TTGCAGAC": "UMI_29", "CAATGTGG": "UMI_30",
	"ACTAGGAG": "UMI_32",
}

func usage() {
	fmt.Print(
		"extractIdtDuplex - Extract Illumina indexes and IDT duplex barcodes from PacBio BAM file and store the result in read tags.\n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var input *string = flag.String("i", "", "Input bam file.")
	var output *string = flag.String("o", "stdout", "Output bam file.")
	var sampleSheet *string = flag.String("s", "", "A sample sheet as a .csv file with the following header \"Sample,i7,i5\" and corresponding data for each sample in the body of the file")
	flag.Parse()
	flag.Usage = usage

	if *input == "" {
		usage()
		log.Fatalln("ERROR: must input a bam file with -i.")
	}

	if *sampleSheet == "" {
		usage()
		log.Fatalln("ERROR: must input a sample sheet with -s")
	}

	var out io.Writer
	switch *output {
	case "stdout":
		out = os.Stdout
	case "stderr":
		out = os.Stderr
	default:
		outfile, err := os.Create(*output)
		exception.PanicOnErr(err)
		defer outfile.Close()
		out = outfile
	}

	extractIdtDuplex(*input, out, *sampleSheet)
}

func extractIdtDuplex(input string, output io.Writer, sampleSheet string) {
	bamChan, bamHeader := sam.GoReadToChan(input)

	sampleIndexes := makeIndexMap(sampleSheet)
	addReadGroupsToHeader(sampleIndexes, &bamHeader)
	out := sam.NewBamWriter(output, bamHeader)

	readsPerSample := make(map[string]int)
	for _, samp := range sampleIndexes {
		readsPerSample[samp] = 0
	}

	var err error
	var isRev bool
	var i7Idx, i5Idx, i7Umi, i5Umi, stringSeq, consensusSampleIdx, sample, sortedUmiPair string
	var numReads, numSuccess int
	var tmpSeq []dna.Base = make([]dna.Base, 0, 100)
	var i7ThreePrimeFlank, i7FivePrimeFlank, i5ThreePrimeFlank, i5FivePrimeFlank string
	for b := range bamChan {
		numReads++
		stringSeq = dna.BasesToString(b.Seq)

		switch {
		case strings.HasPrefix(stringSeq, fwdStart):
			isRev = false
			i7ThreePrimeFlank, i7FivePrimeFlank, i5ThreePrimeFlank, i5FivePrimeFlank = i7ThreePrimeFlankFwd, i7FivePrimeFlankFwd, i5ThreePrimeFlankFwd, i5FivePrimeFlankFwd

		case strings.HasPrefix(stringSeq, revStart):
			isRev = true
			i7ThreePrimeFlank, i7FivePrimeFlank, i5ThreePrimeFlank, i5FivePrimeFlank = i7ThreePrimeFlankRev, i7FivePrimeFlankRev, i5ThreePrimeFlankRev, i5FivePrimeFlankRev

		default: // malformed adapter
			continue
		}

		i7Idx = pullSeqFromFlanks(&b, tmpSeq, isRev, stringSeq, i7FivePrimeFlank, i7ThreePrimeFlank)
		i5Idx = pullSeqFromFlanks(&b, tmpSeq, isRev, stringSeq, i5FivePrimeFlank, i5ThreePrimeFlank)
		if !isRev { // FWD
			i7Umi = consensus(pullSeqBeforeFlank(&b, tmpSeq, !isRev, stringSeq, i7FivePrimeFlankFwd), umiMap)
			i5Umi = consensus(pullSeqAfterFlank(&b, tmpSeq, isRev, stringSeq, i5ThreePrimeFlankFwd), umiMap)
		} else { // REV
			i7Umi = consensus(pullSeqAfterFlank(&b, tmpSeq, !isRev, stringSeq, i7ThreePrimeFlankRev), umiMap)
			i5Umi = consensus(pullSeqBeforeFlank(&b, tmpSeq, isRev, stringSeq, i5FivePrimeFlankRev), umiMap)
		}

		if i7Idx == "" || i5Idx == "" || i7Umi == "" || i5Umi == "" {
			continue
		}

		if !isRev { // FWD
			trimRead(&b, stringSeq, i5ThreePrimeFlank, i7FivePrimeFlank)
		} else { // REV
			trimRead(&b, stringSeq, i7ThreePrimeFlank, i5FivePrimeFlank)
		}

		consensusSampleIdx = consensus(i7Idx+"-"+i5Idx, sampleIndexes)
		sample = sampleIndexes[consensusSampleIdx]
		if _, ok := readsPerSample[sample]; ok {
			readsPerSample[sample]++
		}

		// sortedUmiPair organizes the i5 and i7 UMIs such that sortedUmiPair is the same for reads that
		// have the same two UMIs, but are swapped between i5 and i7 (i.e. strands of a duplex).
		// This is useful for visualizing read families in downstream analysis.
		if i5Umi < i7Umi {
			sortedUmiPair = i5Umi + "-" + i7Umi
		} else {
			sortedUmiPair = i7Umi + "-" + i5Umi
		}

		err = sam.ParseExtra(&b)
		exception.PanicOnErr(err)
		if b.Extra != "" {
			b.Extra += "\t"
		}
		b.Extra += fmt.Sprintf("RG:Z:%s\tAL:Z:%s\tBC:Z:%s\tBF:Z:%s\tBR:Z:%s", sample, sortedUmiPair, i5Umi+"-"+i5Umi, i5Umi, i7Umi)
		sam.WriteToBamFileHandle(out, b, 0)

		numSuccess++
	}

	log.Println("Sample\tReads")
	var lines []string
	for samp, num := range readsPerSample {
		lines = append(lines, fmt.Sprintf("%s\t%d", samp, num))
	}
	slices.Sort(lines)
	log.Println(strings.Join(lines, "\n"))

	log.Printf("Num Records\t%d\n", numReads)
	log.Printf("Num Matched\t%d\n", numSuccess)
	err = out.Close()
	exception.PanicOnErr(err)
}

func pullSeqFromFlanks(b *sam.Sam, tmpSeq []dna.Base, isRev bool, stringSeq, fivePrimeFlank, threePrimeFlank string) string {
	idxStart := strings.Index(stringSeq, fivePrimeFlank)
	idxEnd := strings.Index(stringSeq, threePrimeFlank)
	if idxStart == -1 || idxEnd == -1 {
		return ""
	}
	idxStart += len(fivePrimeFlank)

	if idxEnd-idxStart > 20 || idxEnd-idxStart < 0 {
		return ""
	}

	tmpSeq = tmpSeq[0 : idxEnd-idxStart]
	copy(tmpSeq, b.Seq[idxStart:idxEnd])
	if isRev {
		dna.ReverseComplement(tmpSeq)
	}

	return dna.BasesToString(tmpSeq)
}

func pullSeqAfterFlank(b *sam.Sam, tmpSeq []dna.Base, isRev bool, stringSeq, flank string) string {
	idx := strings.Index(stringSeq, flank)
	if idx == -1 {
		return ""
	}
	idx += len(flank)
	tmpSeq = tmpSeq[0:8]
	copy(tmpSeq, b.Seq[idx:idx+8])
	if isRev {
		dna.ReverseComplement(tmpSeq)
	}
	return dna.BasesToString(tmpSeq)
}

func pullSeqBeforeFlank(b *sam.Sam, tmpSeq []dna.Base, isRev bool, stringSeq, flank string) string {
	idx := strings.Index(stringSeq, flank)
	if idx < 8 {
		return ""
	}
	tmpSeq = tmpSeq[0:8]
	copy(tmpSeq, b.Seq[idx-8:idx])
	if isRev {
		dna.ReverseComplement(tmpSeq)
	}
	return dna.BasesToString(tmpSeq)
}

func consensus(s string, m map[string]string) string {
	if m[s] != "" {
		return s
	}

	for bc, _ := range m {
		if levenshteinString(s, bc) <= 2 {
			return bc
		}
	}

	return "*"
}

func levenshteinString(s1, s2 string) int {
	if s1 == "" || s2 == "" {
		return numbers.Max(len(s1), len(s2))
	}
	s1len := len(s1)
	s2len := len(s2)
	column := make([]int, len(s1)+1)

	for y := 1; y <= s1len; y++ {
		column[y] = y
	}
	for x := 1; x <= s2len; x++ {
		column[0] = x
		lastkey := x - 1
		for y := 1; y <= s1len; y++ {
			oldkey := column[y]
			var incr int
			if s1[y-1] != s2[x-1] {
				incr = 1
			}

			column[y] = minimum(column[y]+1, column[y-1]+1, lastkey+incr)
			lastkey = oldkey
		}
	}
	return column[s1len]
}

func minimum(a, b, c int) int {
	if a < b {
		if a < c {
			return a
		}
	} else {
		if b < c {
			return b
		}
	}
	return c
}

func makeIndexMap(file string) map[string]string {
	data, err := os.ReadFile(file)
	exception.PanicOnErr(err)

	words := strings.Split(string(data), "\n")
	if words[0] != "Sample,i7,i5" {
		log.Fatalf("ERROR: Sample sheet is malformed. Header line must be \"Sample,i7,i5\", but got \"%s\"\n", words[0])
	}

	ans := make(map[string]string)
	var col []string
	for i := 1; i < len(words); i++ {
		words[i] = strings.Trim(words[i], "\r\n")
		col = strings.Split(words[i], ",")
		ans[col[1]+"-"+col[2]] = col[0]
	}

	return ans
}

func addReadGroupsToHeader(indexMap map[string]string, header *sam.Header) {
	for idx, samp := range indexMap {
		header.Text = append(header.Text, fmt.Sprintf("@RG\tID:%s\tBC:%s", samp, idx))
	}
}

func trimRead(b *sam.Sam, stringSeq, fivePrimeFlank, threePrimeFlank string) {
	var insertStart, insertEnd int
	insertStart = strings.Index(stringSeq, fivePrimeFlank) + 8
	insertEnd = strings.Index(stringSeq, threePrimeFlank) - 8
	b.Seq = b.Seq[insertStart:insertEnd]
	b.Qual = b.Qual[insertStart:insertEnd]
}
