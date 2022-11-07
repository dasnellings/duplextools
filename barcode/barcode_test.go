package barcode

import (
	"fmt"
	"testing"
)

func TestEndTrimIndex(t *testing.T) {
	seq := "ATCGTGGGGAGCCCTGGGAGAAGGCTGCCAACCCACCTTCCTGCTGGACGGTCCCGTTAATGAGCTGCCGCGAGTTCCATTCGCGGCCAGATCTTGTCTCCCCAGCAGCCCACTTTCTGT"
	fmt.Println(seq[getEndTrimIndex(seq):])
}
