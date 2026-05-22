package families

import (
	"fmt"
	"github.com/dasnellings/duplextools/barcode"
	"github.com/vertgenlab/gonomics/sam"
	"golang.org/x/exp/maps"
	"log"
)

func GoAnnotate(reads <-chan sam.Sam, startTolerance int, posMatching, strictPosMatching bool) <-chan sam.Sam {
	out := make(chan sam.Sam, 1000)
	go annotate(reads, out, startTolerance, posMatching, strictPosMatching)
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

func annotate(in <-chan sam.Sam, out chan<- sam.Sam, startTolerance int, posMatching, strictPosMatching bool) {
	m := make(map[string]*family)
	readNameMap := make(map[string]uint)
	var currFamilyId uint
	var id string
	var pairMatched bool
	var familyDetermination uint
	var currFam, prevFam, fam *family
	var bf, br string
	currFam = new(family)
	prevFam = new(family)
	var lastPairStart uint32
	var lastPairChr string
	for r := range in {
		bf, br = barcode.Get(r)
		id = getId(bf, br)
		if id == "" {
			out <- r
			continue
		}

		// clean up any stragglers in map if past last start
		if r.RName != lastPairChr || r.Pos > lastPairStart {
			maps.Clear(readNameMap)
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
			currFamilyId++
			currFam.familyId = currFamilyId
			currFam.chr = r.GetChrom()
			currFam.start, currFam.end = getFamilyBoundaries(&r)
			currFam.mateStart = int(r.PNext) - 1
			m[id] = currFam
		}

		// arbitrarily choose first fwd BC as watson strand
		if currFam.watsonStrandId == "" {
			currFam.watsonStrandId = bf
		}

		// add strand tag
		addStrandTag(&r, bf == currFam.watsonStrandId)

		familyDetermination = findFamily(&r, readNameMap, startTolerance, posMatching, strictPosMatching, currFam, prevFam, &currFamilyId)
		addFamilyTag(&r, familyDetermination)

		if !pairMatched && (r.RNext == r.RName || r.RNext == "=") && r.PNext < r.Pos+5000 && r.PNext > r.Pos-5000 { // only track read names if pair is within 5kb to avoid map getting too large
			readNameMap[r.QName] = familyDetermination
			if r.PNext > lastPairStart || r.RName != lastPairChr {
				lastPairStart = r.PNext
				lastPairChr = r.RName
			}
		}

		if len(readNameMap) > 1000000 {
			log.Println("WARNING: read name map is over 1000000 entries long. Performance will degrade. Please contact bugs@dansnellings.com with information about this run.")
		}

		out <- r
	}
	close(out)
}

func findFamily(r *sam.Sam, readNameMap map[string]uint, startTolerance int, posMatching, strictPosMatching bool, currFam, prevFam *family, currFamilyId *uint) uint {
	// gather booleans used for switch-case below
	familyDetermination, pairMatched := readNameMap[r.QName]

	// check if mate is already assigned to a family
	if pairMatched {
		delete(readNameMap, r.QName)
		return familyDetermination
	}

	switch {
	// check match for current family
	case posMatching && r.RName == currFam.chr && familyBoundariesMatch(r, currFam, 0): // start/end match, probably part of existing family
		familyDetermination = currFam.familyId
		if int(r.PNext)-1 != currFam.mateStart {
			currFam.altMateStarts = append(currFam.altMateStarts, int(r.PNext)-1)
		}

	// check match for previous family
	case posMatching && r.RName == prevFam.chr && familyBoundariesMatch(r, prevFam, 0): // start/end match, probably part of existing family
		familyDetermination = prevFam.familyId
		if int(r.PNext)-1 != prevFam.mateStart {
			prevFam.altMateStarts = append(prevFam.altMateStarts, int(r.PNext)-1)
		}

	// barcode match, part of current family
	case !strictPosMatching && r.RName == currFam.chr && familyBoundariesMatch(r, currFam, startTolerance):
		familyDetermination = currFam.familyId

	// barcode match, part of previous family
	case !strictPosMatching && r.RName == prevFam.chr && familyBoundariesMatch(r, prevFam, startTolerance):
		familyDetermination = prevFam.familyId

	default: // must overwrite existing family
		*currFamilyId++
		prevFam, currFam = currFam, prevFam
		currFam.chr = r.RName
		currFam.start, currFam.end = getFamilyBoundaries(r)
		currFam.mateStart = int(r.PNext) - 1
		currFam.altMateStarts = currFam.altMateStarts[:0] // trim
		currFam.familyId = *currFamilyId
		familyDetermination = currFam.familyId
	}

	return familyDetermination
}

func altStartsMatch(start int, altStarts []int) bool {
	for i := range altStarts {
		if start == altStarts[i] {
			return true
		}
	}
	return false
}

func familyBoundariesMatch(r *sam.Sam, fam *family, tolerance int) bool {
	templateStart, templateEnd := getFamilyBoundaries(r)

	if templateStart == fam.start && templateEnd == fam.end {
		return true
	}

	if templateStart < fam.start-tolerance || templateStart > fam.start+tolerance {
		return false
	}

	if templateEnd < fam.end-tolerance || templateEnd > fam.end+tolerance {
		return false
	}

	return true
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

func getFamilyBoundaries(s *sam.Sam) (start, end int) {
	if sam.IsPosStrand(*s) {
		return s.GetChromStart(), s.GetChromStart() + int(s.TLen)
	} else { // is reverse
		return s.GetChromEnd() + int(s.TLen), s.GetChromEnd()
	}
}
