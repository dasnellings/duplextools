package main

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/sam"
	"testing"
)

func TestCigarClipping(t *testing.T) {
	var s sam.Sam
	s.Cigar = cigar.FromString("100M")
	s.Pos = 50
	clipReadEnds(&s, 3)
	if cigar.ToString(s.Cigar) != "3S94M3S" || s.Pos != 53 {
		t.Error("problem with basic cigar clipping", s.Pos, cigar.ToString(s.Cigar))
	}

	s.Cigar = cigar.FromString("3S94M3S")
	s.Pos = 50
	clipReadEnds(&s, 3)
	if cigar.ToString(s.Cigar) != "6S88M6S" || s.Pos != 53 {
		t.Error("problem with basic cigar clipping", s.Pos, cigar.ToString(s.Cigar))
	}

	s.Cigar = cigar.FromString("3S1I100M1I3S")
	s.Pos = 50
	clipReadEnds(&s, 3)
	if cigar.ToString(s.Cigar) != "6S96M6S" || s.Pos != 52 {
		t.Error("problem with basic cigar clipping", s.Pos, cigar.ToString(s.Cigar))
	}

	s.Cigar = cigar.FromString("3S1I100D100M1I3S")
	s.Pos = 50
	clipReadEnds(&s, 3)
	if cigar.ToString(s.Cigar) != "6S96M6S" || s.Pos != 152 {
		t.Error("problem with basic cigar clipping", s.Pos, cigar.ToString(s.Cigar))
	}

	s.Cigar = cigar.FromString("1M1I1D10M")
	s.Pos = 50
	clipReadEnds(&s, 3)
	if cigar.ToString(s.Cigar) != "3S6M3S" || s.Pos != 53 {
		t.Error("problem with basic cigar clipping", s.Pos, cigar.ToString(s.Cigar))
	}

	s.Cigar = cigar.FromString("10S1M10S")
	s.Pos = 50
	clipReadEnds(&s, 3)
	if cigar.ToString(s.Cigar) != "21S" || s.Pos != 51 {
		t.Error("problem with basic cigar clipping", s.Pos, cigar.ToString(s.Cigar))
	}
}
