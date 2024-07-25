package main

import (
	"bytes"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"net/http"
	"strconv"
	"strings"
)

type exon struct {
	chr       string
	start     int
	end       int
	posStrand bool
	geneId    string
	gene      string
}

const exonsToTargetFile string = "/Users/danielsnellings/Desktop/nf_tapestri_panel_selection/RNA/exonsToTarget.gtf"
const refFile string = "/Users/danielsnellings/resources/hg38.fa"
const flankingBases int = 400

const colorRed = "\033[0;31m"
const colorNone = "\033[0m"

func main() {
	var err error
	ref := fasta.NewSeeker(refFile, "")
	defer ref.Close()
	exonsToTarget := fileio.Read(exonsToTargetFile)
	exons := make([]exon, len(exonsToTarget))
	for i := range exonsToTarget {
		exons[i] = parseExon(exonsToTarget[i])
	}

	geneMap := make(map[string][]exon)
	for i := range exons {
		geneMap[exons[i].geneId] = append(geneMap[exons[i].geneId], exons[i])
	}

	var seq1, seq2, leftIdxSeq, rightIdxSeq, leftSeq, rightSeq string
	var leftJunctionBase, rightJunctionBase string
	var cdnaSeq string
	var cdnaQuery string
	cdnaResponseString := new(bytes.Buffer)
	var transcripts []string
	var cdnaResponse *http.Response
	var leftIdx, rightIdx int
	for geneId, exons := range geneMap {
		if len(exons) != 2 {
			continue
		}

		seq1 = getSeq(exons[0], ref)
		seq2 = getSeq(exons[1], ref)

		switch {
		case !exons[0].posStrand && exons[0].start < exons[1].start:
			leftIdxSeq, rightIdxSeq = seq2[max(len(seq2)-20, 0):], seq1[:min(len(seq1), 20)]
		case !exons[0].posStrand && exons[0].start > exons[1].start:
			leftIdxSeq, rightIdxSeq = seq1[max(len(seq1)-20, 0):], seq2[:min(len(seq2), 20)]
		case exons[0].posStrand && exons[0].start < exons[1].start:
			leftIdxSeq, rightIdxSeq = seq1[max(len(seq1)-20, 0):], seq2[:min(len(seq2), 20)]
		case exons[0].posStrand && exons[0].start > exons[1].start:
			leftIdxSeq, rightIdxSeq = seq2[max(len(seq2)-20, 0):], seq1[:min(len(seq1), 20)]
		}

		cdnaQuery = fmt.Sprintf("https://rest.ensembl.org/sequence/id/%s?content-type=text/x-fasta;type=cdna;multiple_sequences=1", strings.Split(geneId, ".")[0])
		cdnaResponse, err = http.Get(cdnaQuery)
		exception.PanicOnErr(err)
		cdnaResponseString.Reset()
		_, err = cdnaResponseString.ReadFrom(cdnaResponse.Body)
		exception.PanicOnErr(err)
		transcripts = readFasta(strings.Split(cdnaResponseString.String(), "\n"))

		cdnaSeq = ""
		for i := range transcripts {
			leftIdx = strings.LastIndex(transcripts[i], leftIdxSeq)
			rightIdx = strings.LastIndex(transcripts[i], rightIdxSeq)
			if leftIdx != -1 && rightIdx != -1 {
				cdnaSeq = transcripts[i]
				break
			}
		}

		leftIdx += 19 // because it indexes the start of the seq. Adding 19 gives us the base on the 5' end of the junction
		leftSeq = cdnaSeq[max(leftIdx-flankingBases, 0):leftIdx]
		leftJunctionBase = string(cdnaSeq[leftIdx])
		rightJunctionBase = string(cdnaSeq[rightIdx])
		rightSeq = cdnaSeq[rightIdx+1 : min(rightIdx+flankingBases+1, len(cdnaSeq))]

		// geneName, geneId, leftSeq[5'base 3'base]rightSeq
		fmt.Printf("%s\t%s\t%s[%s%s]%s\n", exons[0].gene, exons[0].geneId, leftSeq, leftJunctionBase, rightJunctionBase, rightSeq)
	}
}

func getSeq(e exon, ref *fasta.Seeker) string {
	seq, err := fasta.SeekByName(ref, e.chr, e.start, e.end)
	exception.PanicOnErr(err)
	if !e.posStrand {
		dna.ReverseComplement(seq)
	}
	dna.AllToUpper(seq)
	return dna.BasesToString(seq)
}

// chrX	HAVANA	exon	41695555	41697277	.	+	.	gene_id "ENSG00000171659.13"; transcript_id "ENSG00000171659.13"; gene_type "protein_coding"; gene_name "GPR34"; transcript_type "protein_coding"; transcript_name "GPR34"; level 2; havana_gene "OTTHUMG00000021375.2"; exon_id "ENSG00000171659.13_3; exon_number 3";
func parseExon(s string) exon {
	var ans exon
	var err error
	words := strings.Split(s, "\t")
	ans.chr = words[0]
	ans.start, err = strconv.Atoi(words[3])
	exception.PanicOnErr(err)
	ans.end, err = strconv.Atoi(words[4])
	exception.PanicOnErr(err)
	ans.start-- // convert to 0-base left closed, right open
	switch words[6] {
	case "+":
		ans.posStrand = true
	case "-":
		ans.posStrand = false
	default:
		log.Panicf("unrecognized strand: '%s'", words[6])
	}

	ans.geneId = strings.Split(words[8], "\"")[1]
	ans.gene = strings.Split(strings.Split(words[8], "gene_name")[1], "\"")[1]
	return ans
}

func readFasta(lines []string) []string {
	var ans []string
	buf := new(bytes.Buffer)

	var line string
	for _, line = range lines {
		if strings.HasPrefix(line, ">") {
			if buf.Len() > 0 {
				ans = append(ans, buf.String())
				buf.Reset()
			}
			continue
		}

		buf.WriteString(line)
	}

	if buf.Len() > 0 {
		ans = append(ans, buf.String())
	}

	return ans
}

func min(i, j int) int {
	if i < j {
		return i
	}
	return j
}

func max(i, j int) int {
	if i > j {
		return i
	}
	return j
}
