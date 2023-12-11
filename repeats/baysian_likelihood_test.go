package repeats

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"path"
	"strconv"
	"strings"
	"testing"
)

//func TestPermute(t *testing.T) {
//	s := []int{1, 2, 3, 4, 5}
//	ans := permute(s)
//	for i := range ans {
//		fmt.Println(ans[i])
//	}
//}

//func TestUnique(t *testing.T) {
//	s := []int{1, 1, 7, 2, 3, 3, 7, 3, 4, 6}
//	fmt.Println(unique(s))
//}

var testFile string = "/Users/danielsnellings/Important/Lab/Walsh/scSTR_genotyping/UMD1255_208_cells/data/UMD1255.vcf.gz"
var testArr []int = []int{12, 14, 14, 14, 14, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 22, 22, 23, 24, 24, 24, 24, 26, 26, 26, 26, 26, 26, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 31, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32}

type readVect struct {
	v     vcf.Vcf
	names []string
	reads [][]int
}

type data struct {
	splitReads    [][][]int
	baseGenotype  []int
	bestGenotypes [][]int
	ksValues      [][]float64
}

func TestBestSingleGenotype(t *testing.T) {
	repeatUnitLen := 2
	lambda := 0.8
	stutterProb := 0.1
	stutterMismatchProb := 1e-5
	minReadsPerAllele := 50
	var d data
	var i int
	c := goParseVcfToReadVect(testFile)
	for v := range c {
		if v.v.Chr != "chr12" || v.v.Pos != 10591800 {
			continue
		}
		d = getData(v, repeatUnitLen, lambda, stutterProb, stutterMismatchProb, minReadsPerAllele)
		for i := range v.names {
			fmt.Printf("Sample: %s\tbestGenotype: %v\tAllele Size: %d\tReads: %d\tD: %.2f\n", path.Base(v.names[i]), d.bestGenotypes[i], d.baseGenotype[0], len(d.splitReads[i][0]), d.ksValues[i][0])
			fmt.Printf("Sample: %s\tbestGenotype: %v\tAllele Size: %d\tReads: %d\tD: %.2f\n", path.Base(v.names[i]), d.bestGenotypes[i], d.baseGenotype[1], len(d.splitReads[i][1]), d.ksValues[i][1])
		}

		for i = range v.names {
			if genotypesMatch(d.bestGenotypes[i], d.baseGenotype) {
				continue
			}
		}
	}
}

func getData(r readVect, repeatUnitLen int, lambda, stutterProb, stutterMismatchProb float64, minReadsPerAllele int) data {
	var ans data
	ans.bestGenotypes, ans.baseGenotype, _ = BestSharedGenotype(r.reads, repeatUnitLen, lambda, stutterProb, stutterMismatchProb)
	ans.splitReads, ans.ksValues = TestGenotypeFit(r.reads, ans.baseGenotype, repeatUnitLen, lambda, stutterProb, stutterMismatchProb, minReadsPerAllele)
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
		ans.names[val] = key
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

func simulateData() {
	//p := distuv.Poisson{Lambda: 0.2}

}
