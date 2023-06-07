package repeats

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

func TestBestMatch(t *testing.T) {
	seq := dna.StringToBases("ATAAAAAAAAAGACGACGACACGATG")
	pattern := dna.StringToBases("ACG")
	fmt.Println(bestMatch(seq, pattern))
}
