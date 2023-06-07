package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/duplexTools/context"
	"github.com/dasnellings/duplexTools/strand"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"singleStrandSummary - Get summary information about single-strand variants.\n" +
			"Usage:\n" +
			"singleStrandSummary [options] -i input.vcf -r reference.fasta -g reference.gtf > output.txt\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input VCF file with single-strand variant calls.")
	ref := flag.String("r", "", "Reference FASTA file. Must be indexed (.fai).")
	gtf := flag.String("g", "", "Reference GTF file.")
	pad := flag.Int("pad", 1, "Number of bases to use on either side of variant for context.")
	output := flag.String("o", "stdout", "Output TXT file.")
	flag.Parse()

	if *input == "" || *ref == "" || *gtf == "" {
		usage()
		log.Fatalln("ERROR: must have inputs for -i, -r, and -g")
	}

	singleStrandSummary(*input, *ref, *gtf, *output, *pad)
}

type variantType byte

const (
	snv variantType = iota
	insertion
	deletion
)

func singleStrandSummary(vcfFile, refFile, gtfFile, outFile string, pad int) {
	transcribedStrand := make([]vcf.Vcf, 0, 100)
	untranscribedStrand := make([]vcf.Vcf, 0, 100)
	intergenic := make([]vcf.Vcf, 0, 100)
	var transcribedIns, transcribedDel, untranscribedIns, untranscribedDel, intergenicIns, intergenicDel int

	vcfChan, _ := vcf.GoReadToChan(vcfFile)
	gtfTree := gtf.GenesToIntervalTree(gtf.Read(gtfFile))
	ref := fasta.NewSeeker(refFile, "")
	defer cleanup(ref)

	var overlaps []interval.Interval
	var overlappingGene *gtf.Gene
	for v := range vcfChan {
		if !strings.Contains(v.Info, "SS;Strand=") {
			log.Println("skipped (no strand data in info field):", v)
			continue
		}

		overlaps = interval.Query(gtfTree, v, "any")
		if len(overlaps) == 0 {
			intergenic = append(intergenic, v)
			switch getVarType(v) {
			case insertion:
				intergenicIns++
			case deletion:
				intergenicDel++
			}
			continue
		}

		overlappingGene = overlaps[0].(*gtf.Gene) // TODO how to handle multiple gene overlaps
		if isTranscribedStrand(v, overlappingGene) {
			transcribedStrand = append(transcribedStrand, v)
			switch getVarType(v) {
			case insertion:
				transcribedIns++
			case deletion:
				transcribedDel++
			}
		} else {
			untranscribedStrand = append(untranscribedStrand, v)
			switch getVarType(v) {
			case insertion:
				untranscribedIns++
			case deletion:
				untranscribedDel++
			}
		}
	}

	transcribedStrandContexts := context.GetContextMap(goSendVcf(transcribedStrand), ref, pad, false, true)
	untranscribedStrandContexts := context.GetContextMap(goSendVcf(untranscribedStrand), ref, pad, false, true)
	intergenicContexts := context.GetContextMap(goSendVcf(intergenic), ref, pad, false, true)

	var lines []string
	lines = append(lines, fmt.Sprintf("Sample\tMutation\tContext\tTranscribedStrand\tUntranscribedStrand\tIntergenic"))
	for key := range intergenicContexts {
		for context := range intergenicContexts[key] {
			lines = append(lines, fmt.Sprintf("%s\t%s\t%s\t%d\t%d\t%d", vcfFile, key, context, transcribedStrandContexts[key][context], untranscribedStrandContexts[key][context], intergenicContexts[key][context]))
		}
	}
	slices.Sort(lines[1:])
	lines = append(lines, fmt.Sprintf("%s\t%s\t%s\t%d\t%d\t%d", vcfFile, "INS", "NA", transcribedIns, untranscribedIns, intergenicIns))
	lines = append(lines, fmt.Sprintf("%s\t%s\t%s\t%d\t%d\t%d", vcfFile, "DEL", "NA", transcribedDel, untranscribedDel, intergenicDel))

	out := fileio.EasyCreate(outFile)
	defer cleanup(out)

	//TODO COUNT INDELS BY STRAND
	fmt.Fprintln(out, strings.Join(lines, "\n"))
}

func isTranscribedStrand(v vcf.Vcf, g *gtf.Gene) bool {
	// if strands match, vcf is on coding (untranscribed strand)
	//log.Println(v)
	//log.Println(g.Transcripts[0])
	//log.Println()
	return strand.IsPosStrand(v) != g.Transcripts[0].Strand
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}

func goSendVcf(v []vcf.Vcf) <-chan vcf.Vcf {
	ans := make(chan vcf.Vcf, 100)
	go func() {
		for i := range v {
			ans <- v[i]
		}
		close(ans)
	}()
	return ans
}

func getVarType(v vcf.Vcf) variantType {
	switch {
	case len(v.Ref) == 1 && len(v.Ref) == len(v.Alt[0]):
		return snv

	case len(v.Ref) > len(v.Alt[0]):
		return deletion

	case len(v.Ref) < len(v.Alt[0]):
		return insertion

	default:
		log.Panic("unreachable")
		return 0
	}
}
