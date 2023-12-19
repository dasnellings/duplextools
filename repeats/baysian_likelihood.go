package repeats

import (
	"github.com/vertgenlab/gonomics/numbers"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
	"math"
	"sort"
)

// TestGenotypeFit splits sampleReads by most probable allele of origin (of the alleles in the baseGenotype slice). The output splitReads is organized as splitReads[i][j][k] where i is the sample index, j is the genotype index, and k is the read index.
// The output ksValues contains the D statistic of the kolmogorov-smirnov test organized as ksValues[i][j] where i is the sample index, and j is the genotype index. P values are calculated by taking the reads from samples with bestGenoytpe matching
// the baseGenotype, then bootstrapping the D statistic.
func TestGenotypeFit(sampleReads [][]int, bestGenotypes [][]int, baseGenotype []int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64, minReadsPerAllele int, maxAlleleImbalance float64, bootstrapIterations, bootstrapReadCount int) (splitReads [][][]int, ksValues [][]float64, pValues [][]float64) {
	ksValues = make([][]float64, len(sampleReads))
	splitReads = make([][][]int, len(sampleReads))
	baseGenotypeKsReadsAndWeights := make([][][]float64, len(baseGenotype))
	var i, j int
	for j = range baseGenotype {
		baseGenotypeKsReadsAndWeights[j] = getKsWeights(baseGenotype[j], repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
	}
	for i = range sampleReads {
		splitReads[i] = splitReadsByAllele(sampleReads[i], baseGenotype, repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
		ksValues[i] = make([]float64, len(baseGenotype))
		for j = range baseGenotype {
			if len(splitReads[i][j]) < minReadsPerAllele || failsAlleleImbalance(splitReads[i], maxAlleleImbalance) {
				ksValues[i][j] = 0
			} else {
				ksValues[i][j] = ksTest(splitReads[i][j], baseGenotypeKsReadsAndWeights[j])
			}
		}
	}

	// BOOTSTRAPPING
	var samplesMatchingBaseGenotype []int
	for i = range bestGenotypes {
		if equalGenotypes(bestGenotypes[i], baseGenotype) {
			samplesMatchingBaseGenotype = append(samplesMatchingBaseGenotype, i)
		}
	}
	bootstrapReads := make([]int, bootstrapReadCount)
	bootstrapDvalues := make([][]float64, len(baseGenotype))
	for i = range bootstrapDvalues {
		bootstrapDvalues[i] = make([]float64, bootstrapIterations)
	}
	var sampleIdx, alleleIdx, readIdx int
	for i = 0; i < bootstrapIterations; i++ {
		for alleleIdx = 0; alleleIdx < len(baseGenotype); alleleIdx++ {
			for j = 0; j < bootstrapReadCount; j++ {
				sampleIdx = samplesMatchingBaseGenotype[numbers.RandIntInRange(0, len(samplesMatchingBaseGenotype))]
				readIdx = numbers.RandIntInRange(0, len(splitReads[sampleIdx][alleleIdx]))
				bootstrapReads[j] = sampleReads[sampleIdx][readIdx]
			}
			sort.Sort(sort.IntSlice(bootstrapReads))
			bootstrapDvalues[alleleIdx][i] = ksTest(bootstrapReads, baseGenotypeKsReadsAndWeights[alleleIdx])
		}
	}
	for i = range bootstrapDvalues {
		sort.Sort(sort.Float64Slice(bootstrapDvalues[i]))
	}

	// USE BOOTSTRAPPING TO GET P VALUES
	pValues = make([][]float64, len(sampleReads))
	for i = range pValues {
		pValues[i] = make([]float64, len(baseGenotype))
		for j = range pValues[i] {
			pValues[i][j] = getPvalue(ksValues[i][j], bootstrapDvalues[j])
		}
	}
	return
}

// BestSharedGenotype determines the best overall genotype to describe read lengths observed from multiple samples.
func BestSharedGenotype(sampleReads [][]int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) (bestGenotypes [][]int, bestBaseGenotype []int, logLikelihoods []float64) {
	var possibleGenotypes [][]int
	bestGenotypes = make([][]int, len(sampleReads))
	var genotype, mergeReads []int
	for i := range sampleReads {
		genotype, _ = BestSingleGenotype(sampleReads[i], repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
		bestGenotypes[i] = genotype
		if len(genotype) != 0 {
			possibleGenotypes = append(possibleGenotypes, genotype)
		}
		mergeReads = append(mergeReads, sampleReads[i]...)
	}
	possibleGenotypes = uniqGenotypes(sortGenotypes(possibleGenotypes))
	var bestLogLik float64 = math.Inf(-1)
	var currLogLik float64
	for i := range possibleGenotypes {
		currLogLik = GenotypeLogLikelihood(mergeReads, possibleGenotypes[i], repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
		if currLogLik > bestLogLik {
			bestLogLik = currLogLik
			bestBaseGenotype = possibleGenotypes[i]
		}
	}
	logLikelihoods = make([]float64, len(sampleReads))
	for i := range sampleReads {
		logLikelihoods[i] = GenotypeLogLikelihood(sampleReads[i], bestBaseGenotype, repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
	}

	return
}

// BestSingleGenotype determines the best diploid genotype given the observed read lengths and prior stutter parameters. The possible genotypes are
// defined as all possible diploid permutations of observed read lengths. Each possible genotype is tested and the genotype with the maximum
// likelihood is returned.
func BestSingleGenotype(reads []int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) (genotype []int, logLik float64) {
	possibleAlleles := unique(reads)
	possibleGenotypes := permute(possibleAlleles)
	var maxGenotype []int
	var maxLogLik, currLogLik float64
	maxLogLik = -math.MaxFloat64
	for i := range possibleGenotypes {
		currLogLik = GenotypeLogLikelihood(reads, possibleGenotypes[i], repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
		if currLogLik > maxLogLik {
			maxLogLik = currLogLik
			maxGenotype = possibleGenotypes[i]
		}
	}

	return maxGenotype, maxLogLik
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
		return stutterMismatchProb / float64(lengthDiff)
	}

	// read length matches expected stutter profile
	repeatCountDiff := lengthDiff / repeatUnitLength
	p := distuv.Poisson{Lambda: poissonLambda}
	return (stutterProb / 2) * p.Prob(float64(repeatCountDiff-1)) // minus 1 because we want to use the full distribution, the length not changing is not modeled in the poisson.
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

func splitReadsByAllele(reads []int, genotype []int, repeatUnitLength int, poissonLambda, stutterProb, stutterMismatchProb float64) [][]int {
	ans := make([][]int, len(genotype))
	var i, j int
	var curr, maxLik float64
	var assignedGenotype int
	var tiedGenotypes []int
	var prevRead int
	for i = range reads {
		if reads[i] != prevRead {
			//sum = 0
			maxLik = math.Inf(-1)
			tiedGenotypes = tiedGenotypes[:0]
			for j = range genotype {
				curr = likelihoodPerAllele(reads[i], genotype[j], repeatUnitLength, poissonLambda, stutterProb, stutterMismatchProb)
				if curr > maxLik {
					maxLik = curr
					tiedGenotypes = tiedGenotypes[:0]
					tiedGenotypes = append(tiedGenotypes, j)
					continue
				}
				if curr == maxLik {
					tiedGenotypes = append(tiedGenotypes, j)
				}
			}
		}
		prevRead = reads[i]
		assignedGenotype = tiedGenotypes[numbers.RandIntInRange(0, len(tiedGenotypes))]
		ans[assignedGenotype] = append(ans[assignedGenotype], reads[i])
	}
	return ans
}

func sortGenotypes(s [][]int) [][]int {
	sort.Slice(s, func(i, j int) bool {
		var k int
		for k = range s[i] {
			if s[i][k] < s[j][k] {
				return true
			}
			if s[i][k] > s[j][k] {
				return false
			}
		}
		return true
	})
	return s
}

func uniqGenotypes(s [][]int) [][]int {
	if len(s) < 2 {
		return s
	}
	ans := make([][]int, 1)
	ans[0] = s[0]
	for i := 1; i < len(s); i++ {
		if !equalGenotypes(s[i], s[i-1]) {
			ans = append(ans, s[i])
		}
	}
	return ans
}

func equalGenotypes(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func ksTest(reads []int, alleleReadsAndWeights [][]float64) float64 {
	var D float64
	if len(reads) > 1 {
		D = stat.KolmogorovSmirnov(intsToFloats(reads), nil, alleleReadsAndWeights[0], alleleReadsAndWeights[1])
	}
	return D
}

// TODO Optimize
func getKsWeights(alleleSize, repeatUnitLen int, lambda, stutterProb, stutterMismatchProb float64) [][]float64 {
	ans := make([][]float64, 2)
	ans[0] = make([]float64, 100)
	ans[1] = make([]float64, 100)

	for i := 0; i < len(ans[0]); i++ {
		ans[0][i] = float64(i)
		ans[1][i] = likelihoodPerAllele(i, alleleSize, repeatUnitLen, lambda, stutterProb, stutterMismatchProb)
	}

	return ans
}

func intsToFloats(s []int) []float64 {
	ans := make([]float64, len(s))
	for i := range s {
		ans[i] = float64(s[i])
	}
	return ans
}

func genotypesMatch(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func failsAlleleImbalance(genotypeReads [][]int, maxAlleleImbalance float64) bool {
	var min, max, curr float64
	min = 1
	var sum int
	for i := range genotypeReads {
		sum += len(genotypeReads[i])
	}
	for i := range genotypeReads {
		curr = float64(len(genotypeReads[i])) / float64(sum)
		if curr < min {
			min = curr
		}
		if curr > max {
			max = curr
		}
	}
	if max > maxAlleleImbalance || min < 1-maxAlleleImbalance {
		return true
	}
	return false
}

func getPvalue(value float64, bootstrapValues []float64) float64 {
	var idx int
	for idx = range bootstrapValues {
		if value <= bootstrapValues[idx] {
			break
		}
	}
	return float64(len(bootstrapValues)-idx) / float64(len(bootstrapValues))
}