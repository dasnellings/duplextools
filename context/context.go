package context

import (
	"github.com/dasnellings/MCS_MS/strand"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func GetContextMap(vcfChan <-chan vcf.Vcf, ref *fasta.Seeker, pad int, mergeComplements bool, considerStrand bool) map[string]map[string]int {
	ans := initMap(pad)
	for v := range vcfChan {
		vcfContext(v, ans, ref, pad, considerStrand)
	}

	if mergeComplements {
		complementMerge(ans)
	}

	return ans
}

func vcfContext(v vcf.Vcf, m map[string]map[string]int, ref *fasta.Seeker, pad int, considerStrand bool) {
	var topKey, botKey string
	var keyFound bool
	var seq []dna.Base
	var numIncluded int
	var err error
	var needsComplement bool

	// exclude multiallelic and indels
	if !vcf.IsBiallelic(v) || !vcf.IsSubstitution(v) || v.Pos == 1 {
		return
	}

	if considerStrand && !strand.IsPosStrand(v) {
		needsComplement = true
	}

	refBase := v.Ref
	altBase := v.Alt[0]

	if needsComplement {
		refBase = complementString(refBase)
		altBase = complementString(altBase)
	}

	topKey = refBase + ">" + altBase
	// exclude if invalid bases in key
	if _, keyFound = m[topKey]; !keyFound {
		return
	}

	if pad > 0 {
		seq, err = fasta.SeekByName(ref, v.Chr, (v.Pos-1)-pad, (v.Pos-1)+pad+1)
		dna.AllToUpper(seq)
		if needsComplement {
			dna.ReverseComplement(seq)
		}
	} else {
		seq = dna.StringToBases(refBase)
	}
	// exclude if ref sequence does not match
	if err != nil || seq[pad] != dna.StringToBase(refBase) {
		return
	}

	botKey = dna.BasesToString(seq)
	// exclude if invalid bases in key
	if _, keyFound = m[topKey][botKey]; !keyFound {
		return
	}

	numIncluded++
	m[topKey][botKey]++
}

// initialize map with all valid permutations of nucleotides for given pad size.
func initMap(pad int) map[string]map[string]int {
	m := make(map[string]map[string]int)

	// fill in all top level keys
	b := "ACGT"
	var i, j int
	for i = 0; i < len(b); i++ {
		for j = 0; j < len(b); j++ {
			if i == j {
				continue
			}
			m[string(b[i])+">"+string(b[j])] = make(map[string]int)
		}
	}

	// fill in 2nd level keys
	var pfs []string // pfs == possible flanking sequences
	var fullContextSeq string
	permute(b, "", pad*2, &pfs) // generate all possible flanking sequences
	for key, val := range m {   // loop over keys added in previous loop
		for i := range pfs {
			fullContextSeq = pfs[i][:pad] + string(key[0]) + pfs[i][pad:] // sandwich reference base (key[0]) between possible flanking sequence to generate all possible seqs for given ref basek
			val[fullContextSeq] = 0                                       // enter context sequence in map
		}
	}
	return m
}

// permute generates all possible permutations of characters present in b of length k and stores them in ans.
func permute(b string, s string, k int, ans *[]string) {
	if k == 0 {
		*ans = append(*ans, s)
		return
	}

	for i := 0; i < len(b); i++ {
		permute(b, s+string(b[i]), k-1, ans)
	}
}

func complementMerge(m map[string]map[string]int) {
	merge(m["C>A"], m["G>T"])
	delete(m, "G>T")
	merge(m["C>G"], m["G>C"])
	delete(m, "G>C")
	merge(m["C>T"], m["G>A"])
	delete(m, "G>A")
	merge(m["T>A"], m["A>T"])
	delete(m, "A>T")
	merge(m["T>C"], m["A>G"])
	delete(m, "A>G")
	merge(m["T>G"], m["A>C"])
	delete(m, "A>C")
}

func merge(m1, m2 map[string]int) {
	var b []dna.Base
	var keyFound bool
	for key := range m1 {
		b = dna.StringToBases(key)
		dna.ReverseComplement(b)
		if _, keyFound = m2[dna.BasesToString(b)]; !keyFound {
			log.Panicf("something went horribly wrong\n could not find key %s in map below\n%v\n", dna.BasesToString(b), m2)
		}
		m1[key] += m2[dna.BasesToString(b)]
	}
}

func complementString(s string) string {
	switch s {
	case "A", "a":
		return "T"
	case "C", "c":
		return "G"
	case "G", "g":
		return "C"
	case "T", "t":
		return "A"
	default:
		log.Panicln("unrecognized base:", s)
		return ""
	}
}
