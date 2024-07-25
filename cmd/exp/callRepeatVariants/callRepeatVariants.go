package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/repeats"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/vcf"
	"gonum.org/v1/gonum/stat"
	"log"
	"path"
	"strconv"
	"strings"
)

func usage() {
	fmt.Print(
		"callRepeatVariants - Call somatic variants in repeat regions from the output of genotypeTargetRepeats.\n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var input *string = flag.String("i", "", "Input vcf file generated with genotypeTargetRepeats.")
	var output *string = flag.String("o", "stdout", "Output vcf file.")
	flag.Parse()
	flag.Usage = usage

	if *input == "" {
		usage()
		log.Fatalln("ERROR: must input a vcf file with -i.")
	}

	callRepeatVariants(*input, *output)
}

type readVect struct {
	v     vcf.Vcf
	names []string
	reads [][]int
}

type Data struct {
	v             vcf.Vcf
	names         []string
	splitReads    [][][]int
	baseGenotype  []int
	bestGenotypes [][]int
	ksValues      [][]float64
	pValues       [][]float64
}

func callRepeatVariants(input, output string) {
	repeatUnitLen := 2
	lambda := 0.8
	stutterProb := 0.1
	stutterMismatchProb := 1e-5
	minReadsPerAllele := 50
	maxAlleleImbalance := 0.8
	minMutationD := 0.5
	maxMutationP := 0.05
	minFractionOfCellsWithBaseGenotype := 0.5
	bootstrapIterations := 10000
	bootstrapReadCounts := 100

	var d Data
	var i, j int
	var recordNum, cellsWithBaseGenotype, genotypedCells int
	c := goParseVcfToReadVect(input)
	for v := range c {
		recordNum++
		if recordNum%50 == 0 {
			log.Println(recordNum)
		}

		repeatUnitLen = getRepeatUnitSize(v.v)
		if repeatUnitLen < 2 {
			continue
		}
		d = getData(v, repeatUnitLen, lambda, stutterProb, stutterMismatchProb, minReadsPerAllele, maxAlleleImbalance, bootstrapIterations, bootstrapReadCounts)

		if len(d.baseGenotype) < 2 || d.baseGenotype[1]-d.baseGenotype[0] <= repeatUnitLen {
			continue
		}

		cellsWithBaseGenotype = 0
		genotypedCells = 0
		for i = range v.names {
			if len(d.bestGenotypes[i]) != 0 {
				genotypedCells++
			}
			if genotypesMatch(d.bestGenotypes[i], d.baseGenotype) {
				cellsWithBaseGenotype++
			}
		}
		if float64(cellsWithBaseGenotype)/float64(genotypedCells) < minFractionOfCellsWithBaseGenotype {
			continue
		}

		for i = range v.names {
			if genotypesMatch(d.bestGenotypes[i], d.baseGenotype) {
				continue
			}
			for j = range d.ksValues[i] {
				if d.ksValues[i][j] < minMutationD {
					continue
				}
				if d.pValues[i][j] > maxMutationP {
					continue
				}
				fmt.Printf("Variant: %d\tSample: %s\tbestGenotype: %v\tbaseGenotype: %v\tAllele Size: %d\tReads: %d\tD: %.2f\tP: %.2g\tStdev: %.2f\n", recordNum, path.Base(v.names[i]), d.bestGenotypes[i], d.baseGenotype, d.baseGenotype[j], len(d.splitReads[i][j]), d.ksValues[i][j], d.pValues[i][j], stat.StdDev(intsToFloats(d.splitReads[i][j]), nil))
			}
		}
		if recordNum == 168 {
			PlotData(d)
		}
	}
}

func getData(r readVect, repeatUnitLen int, lambda, stutterProb, stutterMismatchProb float64, minReadsPerAllele int, maxAlleleImbalance float64, bootstrapIterations, bootstrapReadCounts int) Data {
	var ans Data
	ans.v = r.v
	ans.names = r.names
	ans.bestGenotypes, ans.baseGenotype, _ = repeats.BestSharedGenotype(r.reads, repeatUnitLen, lambda, stutterProb, stutterMismatchProb)
	ans.splitReads, ans.ksValues, ans.pValues = repeats.TestGenotypeFit(r.reads, ans.bestGenotypes, ans.baseGenotype, repeatUnitLen, lambda, stutterProb, stutterMismatchProb, minReadsPerAllele, maxAlleleImbalance, bootstrapIterations, bootstrapReadCounts)
	return ans
}

func goParseVcfToReadVect(file string) <-chan readVect {
	out := make(chan readVect, 1000)
	go parseVcfToReadVect(file, out)
	return out
}

func parseVcfToReadVect(file string, out chan<- readVect) {
	records, header := vcf.GoReadToChan(file)
	for v := range records {
		out <- parseVcf(v, header)
	}
	close(out)
}

func parseVcf(v vcf.Vcf, header vcf.Header) readVect {
	var ans readVect
	ans.v = v
	ans.names = make([]string, len(header.Samples))
	for key, val := range header.Samples {
		ans.names[val] = path.Base(key)
	}
	ans.reads = make([][]int, len(header.Samples))

	readLengthIndex := getReadLengthsIndex(v.Format)
	for i := range ans.reads {
		ans.reads[i] = parseSample(v.Samples[i].FormatData[readLengthIndex])
	}
	return ans
}

func parseSample(s string) []int {
	if s == "." {
		return []int{}
	}
	words := strings.Split(s, ";")
	a1 := strings.Split(words[0], ",")
	a2 := strings.Split(words[1], ",")
	if a1[0] == "" {
		a1 = nil
	}
	if a2[0] == "" {
		a2 = nil
	}
	merge := append(a1, a2...)

	reads := make([]int, 0)
	var size, num, i, j int
	var err error
	for i = range merge {
		words = strings.Split(merge[i], "=")
		size, err = strconv.Atoi(words[0])
		exception.PanicOnErr(err)
		num, err = strconv.Atoi(words[1])
		exception.PanicOnErr(err)
		for j = 0; j < num; j++ {
			reads = append(reads, size)
		}
	}
	return reads
}

func getReadLengthsIndex(format []string) int {
	var idx int = len(format) - 1
	if format[idx] == "RL" {
		return idx
	}

	for idx = range format {
		if format[idx] == "RL" {
			return idx
		}
	}

	log.Panic("could not find RL in vcf line")
	return 0
}

func getRepeatUnitSize(v vcf.Vcf) int {
	words := strings.Split(v.Id, "x")
	if len(words) != 2 {
		return 0
	}
	return len(words[1])
}

func genotypesMatch(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func intsToFloats(s []int) []float64 {
	ans := make([]float64, len(s))
	for i := range s {
		ans[i] = float64(s[i])
	}
	return ans
}
