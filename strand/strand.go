package strand

import (
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

func IsPosStrand(v vcf.Vcf) bool {
	if !strings.Contains(v.Info, "Strand=") {
		log.Panicln("told to consider strand, but Strand not found in info field:", v)
		return false
	}
	strand := strings.Split(v.Info, "Strand=")[1][0]
	if strand != '+' && strand != '-' {
		log.Panicln("malformed vcf:", v)
	}
	return strand == '+'
}
