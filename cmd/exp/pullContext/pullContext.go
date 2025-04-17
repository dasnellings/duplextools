package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"net/http"
	"regexp"
	"strconv"
	"strings"
)

var file string = "/Users/danielsnellings/Desktop/umd1255_calls/selected_variants.vcf"
var genome string = "/Users/danielsnellings/resources/hg38.fa"
var databasePfx string = "https://rest.ensembl.org/sequence/id/"
var databaseSfxCDNA string = "?content-type=text/x-fasta;type=cdna"
var databaseSfxCDS string = "/sequence/id/ENST00000288602?content-type=text/x-fasta;type=cds"

func usage() {
	fmt.Print(
		"pullContext - Pull flanking cDNA sequence from variants in VCF file.\n" +
			"Usage:\n" +
			"pullContext [options] -i input.vcf -r ref.fasta > output.txt\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input VCF file with variant calls.")
	ref := flag.String("r", "", "Reference FASTA file. Must be indexed (.fai).")
	numFlankBases := flag.Int("pad", 21, "Number of up/downstream bases to include in output.")
	output := flag.String("o", "stdout", "Output file.")
	flag.Parse()

	out := fileio.EasyCreate(*output)
	pad := *numFlankBases

	if *input == "" || *ref == "" {
		usage()
		log.Fatalln("ERROR: must have inputs for -i, -r")
	}

	seeker := fasta.NewSeeker(*ref, "")
	vChan, _ := vcf.GoReadToChan(*input)

	var err error

	_, err = fmt.Fprintln(out, "Chr\tPos\tRef\tAlt\tGene\tTranscript_ID\tRNA_Mutation\tgDNA_Context\tcDNA_Context_Ref\tcDNA_Context_Alt")
	exception.PanicOnErr(err)

	var seq []dna.Base
	var geneName, transcriptId, rnaMut, refSeq, altSeq string
	for v := range vChan {
		seq, err = fasta.SeekByName(seeker, v.Chr, v.GetChromStart()-pad, v.GetChromEnd()+pad)
		exception.PanicOnErr(err)
		dna.AllToUpper(seq)
		dna.AllToLower(seq[:pad])
		dna.AllToLower(seq[len(seq)-pad:])
		geneName, transcriptId, rnaMut, refSeq, altSeq = getTranscriptInfo(v.Info, len(v.Ref), len(v.Alt[0]), v.Ref, v.Alt[0], pad)
		if geneName == "" {
			continue
		}
		_, err = fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt[0], geneName, transcriptId, rnaMut, dna.BasesToString(seq), refSeq, altSeq)
		exception.PanicOnErr(err)
	}

	err = seeker.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func getTranscriptInfo(s string, refSize, altSize int, ref, alt string, pad int) (string, string, string, string, string) {
	var err error
	//fmt.Println(s)
	words := strings.Split(s, "|")
	geneName := words[3]
	transcriptId := words[6]
	rnaMut := words[9]
	cdnaPos := getCDNAPos(rnaMut, refSize, altSize)

	responseCDNA, err := http.Get(databasePfx + strings.Split(transcriptId, ".")[0] + databaseSfxCDNA)
	exception.PanicOnErr(err)

	responseCDS, err := http.Get(databasePfx + strings.Split(transcriptId, ".")[0] + databaseSfxCDS)
	exception.PanicOnErr(err)

	b1 := new(bytes.Buffer)
	b2 := new(bytes.Buffer)

	_, err = b1.ReadFrom(responseCDNA.Body)
	exception.PanicOnErr(err)
	err = responseCDNA.Body.Close()
	exception.PanicOnErr(err)

	_, err = b2.ReadFrom(responseCDS.Body)
	exception.PanicOnErr(err)
	err = responseCDS.Body.Close()
	exception.PanicOnErr(err)

	var refSeq, altSeq string
	linesCDNA := strings.Split(b1.String(), "\n")
	if !strings.Contains(linesCDNA[0], transcriptId) {
		return "", "", "", "", ""
		log.Fatal("transcripts do not match")
	}
	seqCDNA := strings.Join(linesCDNA[1:], "")
	if strings.Contains(seqCDNA, ">") {
		log.Fatal("transcript corrupted", seqCDNA)
	}

	linesCDS := strings.Split(b2.String(), "\n")
	if !strings.Contains(linesCDS[0], transcriptId) {
		log.Fatal("transcripts do not match")
	}
	seqCDS := strings.Join(linesCDS[1:], "")
	if strings.Contains(seqCDS, ">") {
		log.Fatal("transcript corrupted", seqCDS)
	}

	startCDS := strings.Index(seqCDNA, seqCDS)
	endCDS := startCDS + len(seqCDS)

	//fmt.Println(len(seqCDNA), startCDS, endCDS, cdnaPos)

	var variantPos int
	if strings.Contains(rnaMut, "*") {
		variantPos = endCDS + cdnaPos
	} else {
		variantPos = startCDS + cdnaPos
	}

	var start, end int
	start = variantPos - pad
	end = variantPos + pad
	if refSize > altSize {
		end += refSize - altSize
	} else if refSize < altSize {
		start += 1
		end += 1
	} else {
		end += 1
	}
	if start < 0 || end > len(seqCDNA) {
		return "", "", "", "", ""
	}
	refSeq = seqCDNA[start:end]

	refSeq = strings.ToUpper(refSeq)
	refSeq = strings.ToLower(refSeq[:pad]) + refSeq[pad:]
	refSeq = refSeq[:len(refSeq)-pad] + strings.ToLower(refSeq[len(refSeq)-pad:])
	//fmt.Println(refSeq)

	altSeq = refSeq[:pad]
	cdnaVar := refSeq[pad : len(refSeq)-pad]
	// ADD MIDDLE BIT
	//fmt.Println(cdnaVar)
	switch {
	case refSize == altSize: // SNV
		if cdnaVar == ref {
			altSeq += alt
		} else {
			altSeq += dna.BasesToString(dna.ReverseComplementAndCopy(dna.StringToBases(alt)))
		}
	case refSize > altSize: // Deletion
		// do nothing
	case refSize < altSize: // Insertion
		altSeq += strings.Split(rnaMut, "dup")[1]
	}
	// DONE WITH MIDDLE BIT
	altSeq += refSeq[len(refSeq)-pad:]

	return geneName, transcriptId, rnaMut, refSeq, altSeq
}

func getCDNAPos(s string, refsize, altsize int) int {
	words := strings.Split(s, "_")
	var str string
	if refsize >= altsize || strings.Contains(s, "dup") {
		str = words[0]
	} else {
		str = words[1]
	}
	ans, err := strconv.Atoi(regexp.MustCompile(`[^0-9,-]+`).ReplaceAllString(str, ""))
	exception.PanicOnErr(err)
	if ans > 0 {
		ans--
	}
	return ans
}
