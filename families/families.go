package families

import (
	"fmt"
	"github.com/dasnellings/MCS_MS/barcode"
	"github.com/vertgenlab/gonomics/sam"
)

func GoAnnotate(reads <-chan sam.Sam, startTolerance int) <-chan sam.Sam {
	out := make(chan sam.Sam, 1000)
	go annotate(reads, out, uint32(startTolerance))
	return out
}

type family struct {
	chr      string
	pos      uint32
	familyId uint
}

func annotate(in <-chan sam.Sam, out chan<- sam.Sam, startTolerance uint32) {
	m := make(map[string]family)
	var currFamilyId uint
	var id string
	var currFam family
	for r := range in {
		id = getId(barcode.Get(r))
		if id == "" {
			out <- r
			continue
		}
		currFam = m[id]
		switch {
		case r.RName == currFam.chr && r.Pos <= currFam.pos+startTolerance: // part of existing family
			addFamilyTag(&r, currFam.familyId)
		default: // must overwrite existing family
			currFamilyId++
			currFam = family{
				chr:      r.RName,
				pos:      r.Pos,
				familyId: currFamilyId,
			}
			m[id] = currFam
			addFamilyTag(&r, currFam.familyId)
		}
		out <- r
	}
	close(out)
}

func getId(a, b string) string {
	if a == "" || b == "" {
		return ""
	}
	if a > b {
		return a + "-" + b
	}
	return b + "-" + a
}

func addFamilyTag(s *sam.Sam, famId uint) {
	sam.ParseExtra(s)
	if s.Extra != "" {
		s.Extra += "\t"
	}
	s.Extra += fmt.Sprintf("RF:Z:%d", famId)
}
