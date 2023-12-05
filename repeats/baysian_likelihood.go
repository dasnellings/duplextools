package repeats

import (
	"gonum.org/v1/gonum/stat/distuv"
	"math"
	"sort"
)

// BestGenotypeCombo determines the best overall genotype to describe read lengths observed from multiple samples.
func BestGenotypeCombo(sampleReads [][]int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) {

}

// BestSingleGenotype determines the best diploid genotype given the observed read lengths and prior stutter parameters. The possible genotypes are
// defined as all possible diploid permutations of observed read lengths. Each possible genotype is tested and the genotype with the maximum
// likelihood is returned.
func BestSingleGenotype(reads []int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) (genotype []int, logLik float64) {
	possibleAlleles := unique(reads)
	possibleGenotypes := permute(possibleAlleles)

	var maxGenotypeIdx int
	var maxLogLik, currLogLik float64
	maxLogLik = -math.MaxFloat64
	for i := range possibleGenotypes {
		currLogLik = GenotypeLogLikelihood(reads, possibleGenotypes[i], repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
		if currLogLik > maxLogLik {
			maxLogLik = currLogLik
			maxGenotypeIdx = i
		}
	}

	return possibleGenotypes[maxGenotypeIdx], maxLogLik
}

// GenotypeLogLikelihood returns the log likelihood of the input genotype given the observed read lengths and the prior stutter probability.
func GenotypeLogLikelihood(reads []int, genotype []int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) float64 {
	var logLik float64
	for i := range reads {
		logLik += math.Log(ReadLikelihood(reads[i], genotype, repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb))
	}
	return logLik
}

// ReadLikelihood calculates the likelihood that a read of length readLen was derived from the given genotype and stutter parameters. The returned likelihood
// is the max of the likelihoods that the read derived from any alleles in the genotype slice. This heuristic works best for large separations between alleles.
func ReadLikelihood(readLen int, genotype []int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) float64 {
	var maxLikelihood, currLikelihood float64
	for i := range genotype {
		currLikelihood = likelihoodPerAllele(readLen, genotype[i], repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
		if currLikelihood > maxLikelihood {
			maxLikelihood = currLikelihood
		}
	}
	return maxLikelihood
}

func likelihoodPerAllele(readLen int, alleleSize int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) float64 {
	lengthDiff := abs(readLen - alleleSize)

	// read matches allele size
	if lengthDiff == 0 {
		return 1 - stutterProb
	}

	// read length is not a factor of unit length
	if lengthDiff%repeatUnitLength != 0 {
		return stutterMismatchProb
	}

	// read length matches expected stutter profile
	repeatCountDiff := lengthDiff / repeatUnitLength
	p := distuv.Poisson{Lambda: poissonLambda}
	return (stutterProb / 2) * p.Prob(float64(repeatCountDiff))
}

func abs(i int) int {
	if i < 0 {
		i *= -1
	}
	return i
}

func unique(s []int) []int {
	sort.IntSlice(s).Sort()
	var uniqS []int
	var prev int = -1
	for i := range s {
		if s[i] != prev {
			uniqS = append(uniqS, s[i])
			prev = s[i]
		}
	}
	return uniqS
}

func permute(s []int) [][]int {
	var ans [][]int
	if len(s) < 2 {
		return nil
	}
	for i := 1; i < len(s); i++ {
		ans = append(ans, []int{s[0], s[i]})
	}
	ans = append(ans, permute(s[1:])...)
	return ans
}
