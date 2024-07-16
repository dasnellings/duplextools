package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"mutationMotif - Determine multi-base motif of SNV in VCF file.\n" +
			"Usage:\n" +
			"mutationMotif [options] -i input.vcf -r ref.fasta > output.txt\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input VCF file with variant calls.")
	ref := flag.String("r", "", "Reference FASTA file. Must be indexed (.fai).")
	pad := flag.Int("pad", 20, "Number of up/downstream bases to include in output.")
	output := flag.String("o", "stdout", "Output file.")
	flag.Parse()

	if *input == "" || *ref == "" {
		usage()
		log.Fatalln("ERROR: must have inputs for -i, -r")
	}

	mutationMotif(*input, *output, *ref, *pad)
}

func mutationMotif(input, output, refFile string, pad int) {
	sampleName := strings.TrimRight(input, ".vcf")
	records, _ := vcf.GoReadToChan(input)
	out := fileio.EasyCreate(output)
	defer cleanup(out)

	ref := fasta.NewSeeker(refFile, "")
	defer cleanup(ref)

	ans := make(map[string][][4]dna.Base)
	ans["C>A"] = make([][4]dna.Base, 1+pad*2)
	ans["C>G"] = make([][4]dna.Base, 1+pad*2)
	ans["C>T"] = make([][4]dna.Base, 1+pad*2)
	ans["T>A"] = make([][4]dna.Base, 1+pad*2)
	ans["T>C"] = make([][4]dna.Base, 1+pad*2)
	ans["T>G"] = make([][4]dna.Base, 1+pad*2)
	var seq []dna.Base
	var err error
	var start, end int
	var needsRevComp bool
	var mutation string
	for v := range records {
		if !vcf.IsBiallelic(v) || !vcf.IsSubstitution(v) {
			continue
		}

		needsRevComp = false
		if v.Ref == "A" || v.Ref == "G" {
			needsRevComp = true
			mutation = dna.BaseToString(dna.ComplementSingleBase(dna.StringToBase(v.Ref))) + ">" + dna.BaseToString(dna.ComplementSingleBase(dna.StringToBase(v.Alt[0])))
		} else {
			mutation = v.Ref + ">" + v.Alt[0]
		}
		start = v.Pos - 1 - pad
		end = v.Pos + pad
		seq, err = fasta.SeekByName(ref, v.Chr, start, end)
		exception.PanicOnErr(err)
		if !defineSeq(seq) {
			continue
		}
		if needsRevComp {
			dna.ReverseComplement(seq)
		}
		for i := range seq {
			ans[mutation][i][dna.ToUpper(seq[i])]++
		}
	}

	var s string
	for key := range ans {
		s = sampleName + "\t" + key + "\tA"
		for pos := range ans[key] {
			s += fmt.Sprintf("\t%d", ans[key][pos][dna.A])
		}
		fmt.Println(s)
		s = sampleName + "\t" + key + "\tC"
		for pos := range ans[key] {
			s += fmt.Sprintf("\t%d", ans[key][pos][dna.C])
		}
		fmt.Println(s)
		s = sampleName + "\t" + key + "\tG"
		for pos := range ans[key] {
			s += fmt.Sprintf("\t%d", ans[key][pos][dna.G])
		}
		fmt.Println(s)
		s = sampleName + "\t" + key + "\tT"
		for pos := range ans[key] {
			s += fmt.Sprintf("\t%d", ans[key][pos][dna.T])
		}
		fmt.Println(s)
	}
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}

func defineSeq(s []dna.Base) bool {
	for i := range s {
		if !dna.DefineBase(s[i]) {
			return false
		}
	}
	return true
}
