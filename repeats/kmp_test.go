package repeats

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

func TestBuildKmpFailure(t *testing.T) {
	pattern := dna.StringToBases("ACAGAC")
	failure := BuildKmpFailure(pattern)
	fmt.Println(failure)
}

func TestFindRepeat(t *testing.T) {
	pattern := dna.StringToBases("AATGATGATGCAGTGACGTGG")
	fmt.Println(FindRepeat(pattern))
}
