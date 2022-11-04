package families

import (
	"fmt"
	"github.com/dasnellings/MCS_MS/barcode"
	"github.com/vertgenlab/gonomics/sam"
)

func GoAnnotate(reads <-chan sam.Sam, startTolerance int, posMatching bool) <-chan sam.Sam {
	out := make(chan sam.Sam, 1000)
	go annotate(reads, out, startTolerance, posMatching)
	return out
}

type family struct {
	chr            string
	start          int
	mateStart      int
	altMateStarts  []int
	end            int
	familyId       uint
	watsonStrandId string
}

func annotate(in <-chan sam.Sam, out chan<- sam.Sam, startTolerance int, posMatching bool) {
	m := make(map[string]*family)
	var currFamilyId uint
	var id string
	var currFam, prevFam, fam *family
	var bf, br string
	currFam = new(family)
	prevFam = new(family)
	for r := range in {
		bf, br = barcode.Get(r)
		id = getId(bf, br)
		if id == "" {
			out <- r
			continue
		}

		// retrieve family else make a new one.
		// store previous family for efficiency.
		fam = m[id]
		if fam != currFam {
			if currFam.familyId != 0 {
				prevFam = currFam
			}
			currFam = fam
		}

		if currFam == nil {
			currFam = new(family)
			m[id] = currFam
		}

		// arbitrarily choose first fwd BC as watson strand
		if currFam.watsonStrandId == "" {
			currFam.watsonStrandId = bf
		}

		// add strand tag
		addStrandTag(&r, bf == currFam.watsonStrandId)

		switch {
		case posMatching && r.RName == currFam.chr && (r.GetChromStart() == currFam.start || r.GetChromEnd() == currFam.end || int(r.PNext) == currFam.start || r.GetChromStart() == currFam.mateStart): // start/end match, probably part of existing family
			addFamilyTag(&r, currFam.familyId)
			if int(r.PNext)-1 != currFam.mateStart {
				currFam.altMateStarts = append(currFam.altMateStarts, int(r.PNext)-1)
			}

		case posMatching && r.RName == prevFam.chr && (r.GetChromStart() == prevFam.start || r.GetChromEnd() == prevFam.end || int(r.PNext) == prevFam.start || r.GetChromStart() == prevFam.mateStart): // start/end match, probably part of existing family
			addFamilyTag(&r, prevFam.familyId)
			if int(r.PNext)-1 != prevFam.mateStart {
				prevFam.altMateStarts = append(prevFam.altMateStarts, int(r.PNext)-1)
			}

		case posMatching && r.RName == currFam.chr && altStartsMatch(r.GetChromStart(), currFam.altMateStarts):
			addFamilyTag(&r, currFam.familyId)

		case posMatching && r.RName == prevFam.chr && altStartsMatch(r.GetChromStart(), prevFam.altMateStarts):
			addFamilyTag(&r, prevFam.familyId)

		case r.RName == currFam.chr && r.GetChromStart() <= currFam.start+startTolerance: // barcode match, part of existing family
			addFamilyTag(&r, currFam.familyId)

		default: // must overwrite existing family
			currFamilyId++
			prevFam = currFam
			currFam.chr = r.RName
			currFam.start = r.GetChromStart()
			currFam.mateStart = int(r.PNext) - 1
			currFam.altMateStarts = currFam.altMateStarts[:0] // trim
			currFam.end = r.GetChromEnd()
			currFam.familyId = currFamilyId
			addFamilyTag(&r, currFam.familyId)
		}
		out <- r
	}
	close(out)
}

func altStartsMatch(start int, altStarts []int) bool {
	for i := range altStarts {
		if start == altStarts[i] {
			return true
		}
	}
	return false
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
	//sam.ParseExtra(s)
	if s.Extra != "" {
		s.Extra += "\t"
	}
	s.Extra += fmt.Sprintf("RF:Z:%d", famId)
}

func addStrandTag(s *sam.Sam, watsonStrand bool) {
	sam.ParseExtra(s)
	if s.Extra != "" {
		s.Extra += "\t"
	}

	if watsonStrand {
		s.Extra += fmt.Sprintf("RS:Z:W") // watson
	} else {
		s.Extra += fmt.Sprintf("RS:Z:C") // crick
	}
}
