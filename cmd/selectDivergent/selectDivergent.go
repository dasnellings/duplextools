package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"strconv"
	"strings"
)

func usage() {
	fmt.Print(
		"selectDivergent - Output VCF records where there is an allele present that is not in a chosen reference sample.\n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var input *string = flag.String("i", "", "Input VCF file.")
	var refSample *string = flag.String("r", "", "Reference sample name. Sample must be present in input VCF file.")
	var output *string = flag.String("o", "stdout", "Output VCF file.")
	var hetsOnly *bool = flag.Bool("hetOnly", true, "Only consider records where the reference sample is heterozygous.")
	var minGBDiff *int = flag.Int("minGbDiff", 3, "Minimum size difference to count as a separate allele as defined by the GB field in FORMAT. Calculated by HipSTR output. Set to zero to disable.")
	var summary *bool = flag.Bool("summary", true, "Print a summary of divergent sites after run.")
	//var clonal *bool = flag.Bool("clonal", false, "Only output variants present in multiple samples.")
	flag.Parse()
	flag.Usage = usage

	if *input == "" || *refSample == "" {
		usage()
		log.Fatalln("ERROR: must input a VCF file with -i")
	}

	selectDivergent(*input, *output, *refSample, *hetsOnly, *minGBDiff, *summary) //, *clonal)
}

func selectDivergent(input, output, refSamp string, hetsOnly bool, lenDiff int, summary bool) { //, clonal bool) {
	records, header := vcf.GoReadToChan(input)

	var totSites, refSites, hetSites int
	divergentCount := make([]int, len(header.Samples))
	totalCount := make([]int, len(header.Samples))

	if _, ok := header.Samples[refSamp]; !ok {
		log.Fatalf("ERROR: reference sample '%s' not found in input VCF\n", refSamp)
	}
	refSampIdx := header.Samples[refSamp]

	out := fileio.EasyCreate(output)
	vcf.NewWriteHeader(out, header)

	var refAlleles []int16
	var refGB []int
	var alreadyWrote bool
	//var variantSamplesCount int
	for v := range records {
		alreadyWrote = false
		totSites++
		refAlleles = v.Samples[refSampIdx].Alleles
		refGB = getGB(v.Format, v.Samples[refSampIdx].FormatData)
		if len(refAlleles) == 0 {
			continue
		}

		refSites++
		if hetsOnly && refAlleles[0] == refAlleles[1] {
			continue
		}
		hetSites++
		//variantSamplesCount = 0

		for i := range v.Samples {
			if len(v.Samples[i].Alleles) > 0 && v.Samples[i].Alleles[0] != -1 {
				totalCount[i]++
			}
			if i == refSampIdx {
				continue
			}
			if hasDivergent(refAlleles, v.Samples[i].Alleles) && divergentGB(refGB, getGB(v.Format, v.Samples[i].FormatData), lenDiff) {
				divergentCount[i]++
				if !alreadyWrote {
					vcf.WriteVcf(out, v)
				}
			}
		}
	}

	err := out.Close()
	exception.PanicOnErr(err)

	if summary {
		fmt.Fprintf(os.Stderr, "Summary for Reference Sample '%s'\n", refSamp)
		fmt.Fprintf(os.Stderr, "Total Sites Considered:\t%d\n", totSites)
		fmt.Fprintf(os.Stderr, "Total Sites With Ref Sample:\t%d\n", refSites)
		fmt.Fprintf(os.Stderr, "Total Sites Ref Sample Het:\t%d\n", hetSites)
		for key, val := range header.Samples {
			fmt.Fprintf(os.Stderr, "%s Tested Sites:\t%d\n", key, totalCount[val])
			fmt.Fprintf(os.Stderr, "%s Divergent Sites:\t%d\n", key, divergentCount[val])
		}
	}
}

func hasDivergent(ref, test []int16) bool {
	var foundMatch bool
	var i, j int
	for i = range test {
		foundMatch = false
		if test[i] == -1 {
			continue
		}
		for j = range ref {
			if test[i] == ref[j] {
				foundMatch = true
				break
			}
		}
		if foundMatch == false {
			return true
		}
	}
	return false
}

func divergentGB(ref, test []int, minLenDiff int) bool {
	if len(ref) == 0 || len(test) == 0 {
		return false
	}

	var foundMatch bool
	var i, j int
	for i = range test {
		foundMatch = false
		for j = range ref {
			if test[i] <= ref[j]+minLenDiff && test[i] >= ref[j]-minLenDiff {
				foundMatch = true
				break
			}
		}
		if foundMatch == false {
			return true
		}
	}
	return false
}

func getGB(format []string, formatData []string) []int {
	var i int
	for i = range format {
		if format[i] == "GB" {
			break
		}
	}

	gb := strings.Split(formatData[i], "|")
	if len(gb) != 2 {
		return nil
	}

	ans := make([]int, 2)
	var err error
	ans[0], err = strconv.Atoi(gb[0])
	exception.PanicOnErr(err)
	ans[1], err = strconv.Atoi(gb[1])
	exception.PanicOnErr(err)
	return ans
}
