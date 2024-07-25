package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"strconv"
	"strings"
)

func main() {
	infile := flag.String("i", "", "Input VCF file")
	outfile := flag.String("o", "stdout", "Output TSV file")
	gb := flag.Bool("gb", false, "Print GB from format instead of genotype.")
	flag.Parse()

	if *infile == "" {
		flag.PrintDefaults()
	}

	genotypeTable(*infile, *outfile, *gb)
}

func genotypeTable(infile, outfile string, gb bool) {
	var err error
	out := fileio.EasyCreate(outfile)
	records, header := vcf.GoReadToChan(infile)

	h := make([]string, 2+len(header.Samples))
	h[0] = "Chromosome"
	h[1] = "Position"

	for key, val := range header.Samples {
		h[2+val] = key
	}

	fmt.Fprintln(out, strings.Join(h, "\t"))

	s := new(strings.Builder)
	var i, j int
	var currGB []int

	for v := range records {
		s.Reset()
		s.WriteString(fmt.Sprintf("%s\t%d", v.Chr, v.Pos))

		for i = range v.Samples {
			s.WriteString("\t")

			if gb {
				currGB = getGB(v.Format, v.Samples[i].FormatData)
				if len(currGB) > 0 {
					s.WriteString(fmt.Sprintf("%d/%d", currGB[0], currGB[1]))
				} else {
					s.WriteByte('.')
				}
				continue
			}

			for j = range v.Samples[i].Alleles {
				if j > 0 {
					if v.Samples[i].Phase[0] {
						s.WriteByte('|')
					} else {
						s.WriteByte('/')
					}
				}
				s.WriteString(fmt.Sprintf("%d", v.Samples[i].Alleles[j]))
			}
		}
		fmt.Fprintln(out, s.String())
	}

	err = out.Close()
	exception.PanicOnErr(err)
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
