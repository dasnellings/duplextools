package repeats

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
)

// BuildKmpFailure calculates the Knuth-Morris-Pratt failure function for input pattern
// based on https://www.personal.kent.edu/~rmuhamma/Algorithms/MyAlgorithms/StringMatch/kuthMP.htm
func BuildKmpFailure(pattern []dna.Base) []int {
	// failure[i] = length of the longest proper prefix of pattern[0:i] which is also a proper suffix of pattern[0:i]
	// i.e., failure[i] = length of the longest "prefix-suffix" of pattern[0:i]
	failure := make([]int, len(pattern))

	// Length of the previous longest prefix-suffix
	length := 0
	i := 1

	for i < len(pattern) {
		if pattern[i] == pattern[length] {
			failure[i] = length + 1
			length++
			i++
		} else {
			if length > 0 {
				// This is tricky. Consider the example AAACAAAA and i = 7.
				// Also, note that we do not increment i here
				length = failure[length-1]
			} else {
				failure[i] = 0
				i++
			}
		}
	}

	return failure
}

func FindRepeat(seq []dna.Base) (numRepeats int, repeatUnit []dna.Base) {
	failure := BuildKmpFailure(seq)
	fmt.Println(failure)
	lastIndexValue := failure[len(failure)-1]
	length := len(seq) - lastIndexValue

	if length == 0 || length%len(seq) != 0 {
		return len(seq), seq
	}

	substring := seq[:length]
	for i := length; i < len(seq); i += length {
		if dna.CompareSeqsIgnoreCase(seq[i:i+length], substring) != 0 {
			return len(seq), seq
		}
	}

	return len(seq) / length, substring
}
