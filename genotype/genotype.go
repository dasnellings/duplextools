package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

const mcsSharedSequence string = "AGATGTGTATAAGAGACAG"

type minimalRead struct {
	forBarcode string
	revBarcode string
	maxRepeat  int
}

type minimalBed struct {
	chr   string
	start int
	end   int
	name  string
}

func getBarcodes(s sam.Sam) (forward, reverse string) {
	var seq string
	var idxEnd int
	var query any
	var found bool
	var err error
	query, found, err = sam.QueryTag(s, "BF")
	exception.PanicOnErr(err)
	if found {
		seq = query.(string)
		idxEnd = strings.Index(seq, mcsSharedSequence)
		if idxEnd > 0 {
			forward = seq[:idxEnd]
		}
	}
	query, found, err = sam.QueryTag(s, "BR")
	exception.PanicOnErr(err)
	if found {
		seq = query.(string)
		idxEnd = strings.Index(seq, mcsSharedSequence)
		if idxEnd > 0 {
			reverse = seq[:idxEnd]
		}
	}
	return
}

func getLongestRepeat(seq []dna.Base, repeatUnit []dna.Base) int {
	var maxRepeats, currRepeats, currRepeatIdx int
	for i := range seq {
		if seq[i] != repeatUnit[currRepeatIdx] {
			if currRepeats > maxRepeats {
				maxRepeats = currRepeats
			}
			currRepeats = 0
			currRepeatIdx = 0
			continue
		}
		currRepeatIdx++
		if currRepeatIdx == len(repeatUnit) {
			currRepeatIdx = 0
			currRepeats++
		}
	}
	if currRepeats > maxRepeats {
		maxRepeats = currRepeats
	}
	return maxRepeats
}

func getMinimalRead(s sam.Sam, mb minimalBed) minimalRead {
	var mr minimalRead
	var repeatUnit string
	words := strings.Split(mb.name, "x")
	repeatUnit = words[1]
	mr.forBarcode, mr.revBarcode = getBarcodes(s)
	mr.maxRepeat = getLongestRepeat(s.Seq, dna.StringToBases(repeatUnit))
	return mr
}

func getMinimalBed(b bed.Bed) minimalBed {
	var mb minimalBed
	mb.chr = b.Chrom
	mb.start = b.ChromStart
	mb.end = b.ChromEnd
	mb.name = b.Name
	return mb
}

func pileup(alignmentFile string, bedTargets string) map[minimalBed][]minimalRead {
	reads, header := sam.GoReadToChan(alignmentFile)
	if header.Metadata.SortOrder[0] != sam.Coordinate {
		log.Fatal("Input file must be coordinate sorted. check header")
	}
	targets := bed.Read(bedTargets)
	intervalTargets := make([]interval.Interval, len(targets))
	for i := range targets {
		intervalTargets[i] = targets[i]
	}
	tree := interval.BuildTree(intervalTargets)
	var queryResult []interval.Interval
	targetMap := make(map[minimalBed][]minimalRead)
	var mb minimalBed
	var mr minimalRead
	for r := range reads {
		queryResult = interval.Query(tree, r, "d")
		if len(queryResult) == 0 {
			continue
		}
		mb = getMinimalBed(queryResult[0].(bed.Bed))
		mr = getMinimalRead(r, mb)
		if mr.forBarcode != "" && mr.revBarcode != "" {
			targetMap[mb] = append(targetMap[mb], mr)
		}
	}
	return targetMap
}

type allele struct {
	id          string
	watsonReads []minimalRead
	crickReads  []minimalRead
	watsonMode  int
	crickMode   int
}

func splitAlleles(reads []minimalRead) []allele {
	var alleles []allele
	var allelesIdx int
	var id string
	var isWatson bool
	for i := range reads {
		if reads[i].forBarcode < reads[i].revBarcode {
			id = reads[i].forBarcode + reads[i].revBarcode
			isWatson = true
		} else {
			id = reads[i].revBarcode + reads[i].forBarcode
			isWatson = false
		}

		allelesIdx = -1
		for j := range alleles {
			if id == alleles[j].id {
				allelesIdx = j
				break
			}
		}

		if allelesIdx == -1 {
			allelesIdx = len(alleles)
			alleles = append(alleles, allele{id: id})
		}

		if isWatson {
			alleles[allelesIdx].watsonReads = append(alleles[allelesIdx].watsonReads, reads[i])
		} else {
			alleles[allelesIdx].crickReads = append(alleles[allelesIdx].crickReads, reads[i])
		}
	}
	return alleles
}

func filterAlleles(alleles []allele, minAllelicDepth, minStrandedDepth int) []allele {
	var passingIdx []int
	for i := range alleles {
		if len(alleles[i].watsonReads)+len(alleles[i].crickReads) < minAllelicDepth {
			continue
		}
		if len(alleles[i].watsonReads) < minStrandedDepth || len(alleles[i].crickReads) < minStrandedDepth {
			continue
		}
		passingIdx = append(passingIdx, i)
	}
	passingAlleles := make([]allele, len(passingIdx))
	for i := range passingIdx {
		passingAlleles[i] = alleles[passingIdx[i]]
	}
	return passingAlleles
}

func findAllelesMode(alleles []allele) {
	for i := range alleles {
		alleles[i].watsonMode = modeReads(alleles[i].watsonReads)
		alleles[i].crickMode = modeReads(alleles[i].crickReads)
	}
}

func modeReads(r []minimalRead) int {
	m := make(map[int]int)
	for i := range r {
		m[r[i].maxRepeat]++
	}
	var maxCount, maxCountValue int = -1, -1
	for k, v := range m {
		if v > maxCount {
			maxCount = v
			maxCountValue = k
			continue
		}
		if v == maxCount && k < maxCountValue { // on ties, prefer the smaller value
			maxCount = v
			maxCountValue = k
			continue
		}
	}
	return maxCountValue
}

func callGenotypes(target minimalBed, reads []minimalRead, minAllelicDepth, minStrandedDepth int) []allele {
	alleles := splitAlleles(reads)
	alleles = filterAlleles(alleles, minAllelicDepth, minStrandedDepth)
	findAllelesMode(alleles)
	return alleles
}

func debugPrintAlleles(target minimalBed, alleles []allele) {
	sb := new(strings.Builder)
	if len(alleles) > 0 {
		sb.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s", target.chr, target.start, target.end, target.name))
		for i := range alleles {
			sb.WriteByte('\t')
			for j := range alleles[i].watsonReads {
				if j > 0 {
					sb.WriteRune(',')
				}
				sb.WriteString(fmt.Sprintf("%d", alleles[i].watsonReads[j].maxRepeat))
			}
			sb.WriteRune(':')
			for j := range alleles[i].crickReads {
				if j > 0 {
					sb.WriteRune(',')
				}
				sb.WriteString(fmt.Sprintf("%d", alleles[i].crickReads[j].maxRepeat))
			}
			sb.WriteRune(':')
			if alleles[i].watsonMode == alleles[i].crickMode {
				sb.WriteString("PASS")
			} else {
				sb.WriteString("FAIL")
			}
		}
		fmt.Println(sb.String())
	}
}

func printHeader(inputFiles []string) {
	sb := new(strings.Builder)
	sb.WriteString("Chromosome\tStart\tEnd\tRepeat")
	for i := range inputFiles {
		sb.WriteString("\t")
		sb.WriteString(strings.TrimSuffix(inputFiles[i], ".bam"))
	}
	fmt.Println(sb.String())
}

func printAlleles(target minimalBed, alleles []allele) {
	sb := new(strings.Builder)
	var allelesWritten bool
	if len(alleles) > 0 {
		sb.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\t", target.chr, target.start, target.end, target.name))
		allelesWritten = false
		for i := range alleles {
			if alleles[i].watsonMode != -1 && alleles[i].watsonMode == alleles[i].crickMode {
				if allelesWritten {
					sb.WriteRune(',')
				}
				sb.WriteString(fmt.Sprintf("%d", alleles[i].watsonMode))
				allelesWritten = true
			}
		}
		if allelesWritten {
			fmt.Println(sb.String())
		}
	}
}

// inputFiles is a custom type that gets filled by flag.Parse()
type inputFiles []string

// String to satisfy flag.Value interface
func (i *inputFiles) String() string {
	return strings.Join(*i, " ")
}

// Set to satisfy flag.Value interface
func (i *inputFiles) Set(value string) error {
	*i = append(*i, value)
	return nil
}

func main() {
	var inputs inputFiles
	flag.Var(&inputs, "i", "Sam or Bam file with alignments, can be declared more than once")
	targetsFile := flag.String("targets", "", "Bed file with target regions")
	minAllelicDepth := flag.Int("a", 4, "Minimum reads per allele for analysis")
	minStrandedDepth := flag.Int("s", 2, "Minimum reads per strand per allele for analysis")
	flag.Parse()

	if len(inputs) == 0 || *targetsFile == "" {
		flag.PrintDefaults()
		log.Fatalln("ERROR: must declare reads and targets file")
	}

	printHeader(inputs)
	bedTargets := bed.Read(*targetsFile)
	targetMap := make(map[minimalBed][][]int) // first slice is file // second slice is alleles // int is mode of repeat length

	for i := range bedTargets {
		targetMap[getMinimalBed(bedTargets[i])] = make([][]int, len(inputs))
	}

	for i := range inputs {
		targets := pileup(inputs[i], *targetsFile)
		var alleles []allele

		for target, reads := range targets {
			alleles = callGenotypes(target, reads, *minAllelicDepth, *minStrandedDepth)
			for j := range alleles {
				if alleles[j].watsonMode == alleles[j].crickMode {
					targetMap[target][i] = append(targetMap[target][i], alleles[j].watsonMode)
				}
			}
		}
	}

	sb := new(strings.Builder)
	var allelesWritten bool
	for k, v := range targetMap {
		allelesWritten = false
		sb.Reset()
		sb.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s", k.chr, k.start, k.end, k.name))
		for i := range v {
			sb.WriteRune('\t')
			if len(v[i]) == 0 {
				sb.WriteRune('.')
			}
			for j := range v[i] {
				if j > 0 {
					sb.WriteRune(',')
				}
				sb.WriteString(fmt.Sprintf("%d", v[i][j]))
				allelesWritten = true
			}
		}
		if allelesWritten {
			fmt.Println(sb.String())
		}
	}
}
