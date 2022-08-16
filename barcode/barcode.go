package barcode

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

const McsSharedSequence string = "AGATGTGTATAAGAGACAG"
const McsSharedSequenceRevComp string = "CTGTCTCTTATACACATCT"
const McsB1 string = "GGCACCGAAAA"
const McsB2 string = "CTCGGCGATAAA"
const McsB3 string = "GGTGGAGCATAA"
const McsB4 string = "CGAGCGCATTAA"
const McsB5 string = "AGCCCGGTTATA"
const McsB6 string = "TCGGCACCAATA"
const McsB7 string = "GCCTGTGGATTA"
const McsB8 string = "GCGACCCTTTTA"
const McsB9 string = "GCATGCGGTAAT"
const McsB10 string = "GCGTTGCCATAT"
const McsB11 string = "GGCCGCATTTAT"
const McsB12 string = "ACCGCCTCTATT"
const McsB13 string = "CCGTGCCAAAAT"
const McsB14 string = "TCTCCGGGAATT"
const McsB15 string = "CCGCGCTTATTT"
const McsB16 string = "CTGAGCTCGTTTT"

var Barcodes map[string]bool = map[string]bool{
	McsB1:  true,
	McsB2:  true,
	McsB3:  true,
	McsB4:  true,
	McsB5:  true,
	McsB6:  true,
	McsB7:  true,
	McsB8:  true,
	McsB9:  true,
	McsB10: true,
	McsB11: true,
	McsB12: true,
	McsB13: true,
	McsB14: true,
	McsB15: true,
	McsB16: true,
}

func Get(s sam.Sam) (forward, reverse string) {
	//var seq string
	//var idxEnd int
	var query any
	var found bool
	var err error
	query, found, err = sam.QueryTag(s, "BF")
	exception.PanicOnErr(err)
	if found {
		forward = query.(string)
		//seq = query.(string)
		//idxEnd = strings.Index(seq, mcsSharedSequence)
		//if idxEnd > 0 {
		//	forward = seq[:idxEnd]
		//}
	}
	query, found, err = sam.QueryTag(s, "BR")
	exception.PanicOnErr(err)
	if found {
		reverse = query.(string)
		//seq = query.(string)
		//idxEnd = strings.Index(seq, mcsSharedSequence)
		//if idxEnd > 0 {
		//	reverse = seq[:idxEnd]
		//}
	}
	return
}

func Trim(fq *fastq.Fastq) {
	s := dna.BasesToString(fq.Seq)
	templateStart := strings.LastIndex(s, McsSharedSequence) + len(McsSharedSequence)
	templateEnd := strings.Index(s, McsSharedSequenceRevComp)
	if templateStart == -1 {
		templateStart = 0
	}
	if templateEnd == -1 {
		templateEnd = len(s)
	}

	if templateStart > templateEnd {
		templateStart = 0
		templateEnd = 0
	}
	fq.Seq = fq.Seq[templateStart:templateEnd]
	fq.Qual = fq.Qual[templateStart:templateEnd]
}

func Extract(seq []dna.Base) string {
	s := dna.BasesToString(seq)
	idx := strings.Index(s, McsSharedSequence)
	if idx == -1 {
		return "*"
	}
	bc := s[:idx]
	return Consensus(bc)
}

// Consensus returns the matching MCS barcode. When an exact match is not found,
// Levenshtein distance is used to predict what the most probable barcode is.
// If Levenshtein distance is >2, then Consensus returns "*".
func Consensus(s string) string {
	if Barcodes[s] {
		return s
	}

	for bc, _ := range Barcodes {
		if levenshteinString(s, bc) <= 2 {
			return bc
		}
	}

	return "*"
}

func revComp(base string) string {
	ans := make([]byte, len(base))
	var j int
	for i := len(base) - 1; i >= 0; i-- {
		switch base[i] {
		case 'A', 'a':
			ans[j] = 'T'
		case 'C', 'c':
			ans[j] = 'G'
		case 'G', 'g':
			ans[j] = 'C'
		case 'T', 't':
			ans[j] = 'A'
		case 'N', 'n':
			ans[j] = 'N'
		default:
			log.Panicf("ERROR: unrecognized base '%s'", base)
		}
		j++
	}
	return string(ans)
}

func levenshtein(s1, s2 []dna.Base) int {
	if len(s1) == 0 || len(s2) == 0 {
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

func minLevenshtein() int {
	barcodes := []string{McsB1,
		McsB2,
		McsB3,
		McsB4,
		McsB5,
		McsB6,
		McsB7,
		McsB8,
		McsB9,
		McsB10,
		McsB11,
		McsB12,
		McsB13,
		McsB14,
		McsB15,
		McsB16}

	var min, curr int = 10, 0
	for i := range barcodes {
		for j := range barcodes {
			if j == i {
				continue
			}
			curr = levenshteinString(barcodes[i], barcodes[j])
			fmt.Printf("Curr:\t%d\t%d\t%d\n", i, j, curr)
			if curr < min {
				min = curr
			}
		}
	}
	fmt.Println("Min:", min)
	return min
}
