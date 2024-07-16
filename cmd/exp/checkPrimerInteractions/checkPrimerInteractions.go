package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"golang.org/x/exp/slices"
	"os"
	"path"
	"sort"
	"strings"
	"sync"
)

var fwdAdapter string = "GTACTCGCAGTAGTC"
var revAdapter string = "GTGATACACGACTATGAGCGCTA"
var optimizaitonIterations int = 8
var superIterations int = 30
var primersFolder string = "/Users/danielsnellings/Important/Lab/Walsh/Neurofibromatosis/nf_tapestri_panel_selection/RNA/primer3_design/primer_files"
var dnaInfile string = "/Users/danielsnellings/Important/Lab/Walsh/Neurofibromatosis/nf_tapestri_panel_selection/Tapestri-Designer-results-7538/Designer/7538-design-summary.csv"
var rnaInfile string = "/Users/danielsnellings/Important/Lab/Walsh/Neurofibromatosis/nf_tapestri_panel_selection/RNA/primer3_design/selected_primers_trim.tsv"
var outfile string = "/Users/danielsnellings/Important/Lab/Walsh/Neurofibromatosis/nf_tapestri_panel_selection/RNA/primer3_design/final_rna_primer_set.txt"
var rnaRevBasesToCheck int = 5
var exactMatchBasesToCheck int = 8
var maxShuffleStart int = 5

type iterationResult struct {
	rnaConflicts        int
	exactMatchConflicts int
	rnaIds              []string
	rnaFwd              []string
	rnaRev              []string
}

func main() {
	//allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev := selectOptimalPrimers()
	allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev := readOutput(outfile)

	allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev = removeProblematicPrimers(allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev)

	out, err := os.Create(outfile)
	exception.PanicOnErr(err)
	defer out.Close()

	_, err = out.WriteString("id\tforward\treverse\n")
	exception.PanicOnErr(err)
	for i := range allIds {
		_, err = out.WriteString(fmt.Sprintf("%s\t%s\t%s\n", allIds[i], allFwd[i], allRev[i]))
		exception.PanicOnErr(err)
	}
}

func selectOptimalPrimers() (allIds, allForPrimers, allRevPrimers, dnaIds, dnaForPrimers, dnaRevPrimers, rnaIds, rnaForPrimers, rnaRevPrimers []string) {
	// read in primers from tapestri designer and primer3
	primersToChoose := readPrimerFiles(primersFolder)
	dnaIds, dnaForPrimers, dnaRevPrimers = readDnaPrimers(dnaInfile)
	rnaIds, rnaForPrimers, rnaRevPrimers = readRnaPrimers(rnaInfile)
	allIds = append(dnaIds, rnaIds...)
	allForPrimers = append(dnaForPrimers, rnaForPrimers...)
	allRevPrimers = append(dnaRevPrimers, rnaRevPrimers...)
	dnaIds = allIds[:len(dnaIds)]
	dnaForPrimers = allForPrimers[:len(dnaForPrimers)]
	dnaRevPrimers = allRevPrimers[:len(dnaRevPrimers)]
	rnaIds = allIds[len(dnaIds):]
	rnaForPrimers = allForPrimers[len(dnaForPrimers):]
	rnaRevPrimers = allRevPrimers[len(dnaRevPrimers):]

	iterationResults := make(chan iterationResult, superIterations)
	wg := new(sync.WaitGroup)
	for i := 0; i < superIterations; i++ {
		// arrange slices so there is allN slice and the rnaN and dnaN slices are subslices
		currAllIds := make([]string, len(allIds))
		copy(currAllIds, allIds)
		currAllForPrimers := make([]string, len(allForPrimers))
		copy(currAllForPrimers, allForPrimers)
		currAllRevPrimers := make([]string, len(allRevPrimers))
		copy(currAllRevPrimers, allRevPrimers)

		currDnaIds := currAllIds[:len(dnaIds)]
		currDnaForPrimers := currAllForPrimers[:len(dnaForPrimers)]
		currDnaRevPrimers := currAllRevPrimers[:len(dnaRevPrimers)]
		currRnaIds := currAllIds[len(dnaIds):]
		currRnaForPrimers := currAllForPrimers[len(dnaForPrimers):]
		currRnaRevPrimers := currAllRevPrimers[len(dnaRevPrimers):]

		if i > 0 {
			shuffleStartingPrimers(currRnaIds, currRnaForPrimers, currRnaRevPrimers, primersToChoose, maxShuffleStart)
		}

		wg.Add(1)
		go runSuperIteration(currAllIds, currAllForPrimers, currAllRevPrimers, currDnaIds, currDnaForPrimers, currDnaRevPrimers, currRnaIds, currRnaForPrimers, currRnaRevPrimers, rnaRevBasesToCheck, exactMatchBasesToCheck, primersToChoose, iterationResults, wg)
	}

	go func() {
		wg.Wait()
		close(iterationResults)
	}()

	var set bool
	var bestRnaConflicts, bestExactMatchConflicts int
	bestRnaForPrimers := make([]string, len(rnaForPrimers))
	bestRnaRevPrimers := make([]string, len(rnaRevPrimers))
	for result := range iterationResults {
		fmt.Printf("Iteration Results:\t%d rna rev conflicts\t%d exact match conflicts\n", result.rnaConflicts, result.exactMatchConflicts)
		if !set || (result.rnaConflicts <= bestRnaConflicts && result.exactMatchConflicts <= bestExactMatchConflicts) {
			set = true
			bestRnaConflicts = result.rnaConflicts
			bestExactMatchConflicts = result.exactMatchConflicts
			copy(bestRnaForPrimers, result.rnaFwd)
			copy(bestRnaRevPrimers, result.rnaRev)
		}
	}

	copy(rnaForPrimers, bestRnaForPrimers)
	copy(rnaRevPrimers, bestRnaRevPrimers)
	fmt.Printf("Best super-iteration values:\n\tRNA reverse match conflicts\t%d\n\tAll exact match conflicts\t%d\n", checkAllRnaConflicts(rnaRevPrimers, rnaRevBasesToCheck), checkExact3PrimeMatch(allIds, allForPrimers, allRevPrimers, exactMatchBasesToCheck))
	return
}

func removeProblematicPrimers(allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev []string) (_allIds, _allFwd, _allRev, _dnaIds, _dnaFwd, _dnaRev, _rnaIds, _rnaFwd, _rnaRev []string) {
	//primersToChoose := readPrimerFiles(primersFolder)
	rnaRevConflicts := checkAllRnaConflicts(rnaRev, rnaRevBasesToCheck)
	exactMatchConflicts := checkExact3PrimeMatch(allIds, allFwd, allRev, exactMatchBasesToCheck)

	var conflictMap map[string]int
	var keys []string
	var key string
	var numRemoved int
	for rnaRevConflicts > 0 {
		conflictMap = getRnaConflictsPerId(rnaIds, rnaRev, rnaRevBasesToCheck)
		keys = make([]string, 0, len(conflictMap))
		for key, _ = range conflictMap {
			keys = append(keys, key)
		}
		sort.Slice(keys, func(i, j int) bool {
			if conflictMap[keys[i]] > conflictMap[keys[j]] {
				return true
			}
			return false
		})

		allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev = removeId(keys[0], allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev)
		numRemoved++
		//optimizeRnaReverse(allIds, allFwd, allRev, rnaIds, rnaFwd, rnaRev, rnaRevBasesToCheck, exactMatchBasesToCheck, primersToChoose)

		rnaRevConflicts = checkAllRnaConflicts(rnaRev, rnaRevBasesToCheck)
		exactMatchConflicts = checkExact3PrimeMatch(allIds, allFwd, allRev, exactMatchBasesToCheck)
		fmt.Printf("Removed %s\t%d rna rev conflicts\t%d exact match conflicts\n", keys[0], rnaRevConflicts, exactMatchConflicts)
	}

	fmt.Printf("Total Genes Removed: %d\n", numRemoved)

	return allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev
}

func runSuperIteration(allIds, allForPrimers, allRevPrimers, dnaIds, dnaForPrimers, dnaRevPrimers, rnaIds, rnaForPrimers, rnaRevPrimers []string, rnaRevBasesToCheck, exactMatchBasesToCheck int, primersToChoose map[string]possiblePrimers, iterationResults chan<- iterationResult, wg *sync.WaitGroup) {
	optimizeRnaReverse(allIds, allForPrimers, allRevPrimers, rnaIds, rnaForPrimers, rnaRevPrimers, rnaRevBasesToCheck, exactMatchBasesToCheck, primersToChoose)
	rnaConflicts := checkAllRnaConflicts(rnaRevPrimers, rnaRevBasesToCheck)
	exactMatchConflicts := checkExact3PrimeMatch(allIds, allForPrimers, allRevPrimers, exactMatchBasesToCheck)
	iterationResults <- iterationResult{
		rnaConflicts:        rnaConflicts,
		exactMatchConflicts: exactMatchConflicts,
		rnaIds:              rnaIds,
		rnaFwd:              rnaForPrimers,
		rnaRev:              rnaRevPrimers,
	}
	wg.Done()
}

// NOTE: names, forPrimers, and revPrimers are a subslice of allIds, allForPrimers, and allRevPrimers respectively.
// Therefore, changes to the first 3 slices, are automatically present in the second 3 slices.
func optimizeRnaReverse(allIds, allForPrimers, allRevPrimers, names, forPrimers, revPrimers []string, rnaRevBasesToCheck, exactMatchBasesToCheck int, primersToChoose map[string]possiblePrimers) {
	var bestFor, bestRev string
	var currConflicts, minConflicts, currFailingExactMatch, minFailingExactMatch int
	var primerOptions possiblePrimers

	fmt.Printf("Optimization iteration 0:\t%d rna rev conflicts\t%d exact match conflicts\n", checkAllRnaConflicts(revPrimers, rnaRevBasesToCheck), checkExact3PrimeMatch(allIds, allForPrimers, allRevPrimers, exactMatchBasesToCheck))

	for k := 0; k < optimizaitonIterations; k++ {
		for i := range names {
			primerOptions = primersToChoose[names[i]]
			for j := range primerOptions.primers {
				forPrimers[i] = fwdAdapter + primerOptions.primers[j][0] // 0 is forward primer, 1 is reverse primer
				revPrimers[i] = revAdapter + primerOptions.primers[j][1]

				currConflicts = 0
				currFailingExactMatch = 0
				for l := range revPrimers {
					currConflicts += checkRnaConflicts(revPrimers[l], revPrimers, rnaRevBasesToCheck)
				}
				currFailingExactMatch = checkExact3PrimeMatch(allIds, allForPrimers, allRevPrimers, exactMatchBasesToCheck)

				if j == 0 || (currConflicts < minConflicts && currFailingExactMatch < minFailingExactMatch+2) || (currConflicts < minConflicts+10 && currFailingExactMatch < minFailingExactMatch) {
					minFailingExactMatch = currFailingExactMatch
					minConflicts = currConflicts
					bestFor = forPrimers[i]
					bestRev = revPrimers[i]
				}
			}
			forPrimers[i] = bestFor
			revPrimers[i] = bestRev
		}
		fmt.Printf("Optimization iteration %d:\t%d rna rev conflicts\t%d exact match conflicts\n", k+1, minConflicts, minFailingExactMatch)
	}

	//fmt.Println(getRnaConflictsPerId(names, revPrimers, rnaRevBasesToCheck))
	//fmt.Println("Problems (checkRNA): ", checkAllRnaConflicts(revPrimers, rnaRevBasesToCheck))
	//fmt.Println("Problems (checkAll): ", checkExact3PrimeMatch(allIds, allForPrimers, allRevPrimers, exactMatchBasesToCheck))
}

func checkRnaConflicts(primer string, primers []string, numBasesToCheck int) int {
	endSeq := primer[len(primer)-numBasesToCheck:]
	compEndSeq := dna.BasesToString(dna.ReverseComplementAndCopy(dna.StringToBases(endSeq)))

	var conflicts int
	for i := range primers {
		if strings.Index(primers[i], compEndSeq) != -1 {
			conflicts++
		}
	}
	return conflicts
}

func checkAllRnaConflicts(primers []string, numBasesToCheck int) int {
	var conflicts int
	for i := range primers {
		conflicts += checkRnaConflicts(primers[i], primers, numBasesToCheck)
	}
	return conflicts
}

func getRnaConflictsPerId(names, primers []string, numBasesToCheck int) map[string]int {
	ans := make(map[string]int)
	var i, j int
	var primer, compEndSeq, endSeq string
	for i, primer = range primers {
		endSeq = primer[len(primer)-numBasesToCheck:]
		compEndSeq = dna.BasesToString(dna.ReverseComplementAndCopy(dna.StringToBases(endSeq)))
		for j = range primers {
			if strings.Index(primers[j], compEndSeq) != -1 {
				ans[names[i]]++
				ans[names[j]]++
			}
		}
	}

	return ans
}

func checkExact3PrimeMatch(names, forPrimers, revPrimers []string, numBasesToCheck int) int {
	primerMap := make(map[string][]string)

	var forward, reverse string
	for i := range names {
		forward = forPrimers[i][len(forPrimers[i])-numBasesToCheck:]
		reverse = revPrimers[i][len(revPrimers[i])-numBasesToCheck:]
		primerMap[forward] = append(primerMap[forward], names[i]+"_for")
		primerMap[reverse] = append(primerMap[reverse], names[i]+"_rev")

		forward = dna.BasesToString(dna.ReverseComplementAndCopy(dna.StringToBases(forward)))
		reverse = dna.BasesToString(dna.ReverseComplementAndCopy(dna.StringToBases(reverse)))
		primerMap[forward] = append(primerMap[forward], names[i]+"_for_revcomp")
		primerMap[reverse] = append(primerMap[reverse], names[i]+"_rev_revcomp")
	}

	var problemKeys int
	for _, val := range primerMap {
		if len(val) > 1 && !allAmpl(val) {
			//fmt.Println(key, len(val), val)
			problemKeys++
		}
	}
	return problemKeys
}

func allAmpl(s []string) bool {
	for i := range s {
		if !strings.HasPrefix(s[i], "AMPL") {
			return false
		}
	}
	return true
}

type possiblePrimers struct {
	geneName string
	primers  [][2]string
}

func readPrimerFiles(dir string) map[string]possiblePrimers {
	files, err := os.ReadDir(dir)
	exception.PanicOnErr(err)

	ans := make(map[string]possiblePrimers)
	var name string
	for i := range files {
		name = strings.Split(files[i].Name(), "_")[0]
		ans[name] = possiblePrimers{geneName: name, primers: readPrimers(path.Join(dir, files[i].Name()))}
	}
	return ans
}

func readPrimers(file string) [][2]string {
	input := fileio.Read(file)
	input = input[1:] // remove header line
	var ans [][2]string
	var words []string
	for i := range input {
		if len(input[i]) == 0 {
			continue
		}
		words = strings.Split(input[i], "\t")
		ans = append(ans, [2]string{words[3], words[4]})
	}
	return ans
}

func readDnaPrimers(file string) (names, forPrimers, revPrimers []string) {
	input := fileio.Read(file)
	input = input[1:] // remove header line

	var words []string
	for i := range input {
		words = strings.Split(input[i], ",")
		names = append(names, words[0])
		forPrimers = append(forPrimers, words[7])
		revPrimers = append(revPrimers, words[8])
	}
	return
}

func readRnaPrimers(file string) (names, forPrimers, revPrimers []string) {
	input := fileio.Read(file)
	input = input[1:] // remove header line

	var words []string
	for _, line := range input {
		words = strings.Split(line, "\t")
		if len(words) < 4 {
			continue // primer3 failed to find primers
		}
		names = append(names, words[0])
		forPrimers = append(forPrimers, fwdAdapter+words[2])
		revPrimers = append(revPrimers, revAdapter+words[3]) // adding in adapter tail for check
	}
	return
}

func shuffleStartingPrimers(ids, fwd, rev []string, primersToChoose map[string]possiblePrimers, maxStart int) {
	var idx int
	var options [][2]string
	for i := range ids {
		options = primersToChoose[ids[i]].primers
		idx = numbers.RandIntInRange(0, min(maxStart, len(options)))
		fwd[i] = options[idx][0]
		rev[i] = options[idx][1]
	}
}

func min(i, j int) int {
	if i < j {
		return i
	}
	return j
}

func readOutput(file string) (allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev []string) {
	input := fileio.Read(file)
	input = input[1:] // remove header line
	var words []string
	for i := range input {
		if len(input[i]) == 0 {
			continue
		}
		words = strings.Split(input[i], "\t")
		if strings.HasPrefix(input[i], "AMPL") {
			dnaIds = append(dnaIds, words[0])
			dnaFwd = append(dnaFwd, words[1])
			dnaRev = append(dnaRev, words[2])
		} else {
			rnaIds = append(rnaIds, words[0])
			rnaFwd = append(rnaFwd, words[1])
			rnaRev = append(rnaRev, words[2])
		}
	}

	allIds = append(dnaIds, rnaIds...)
	allFwd = append(dnaFwd, rnaFwd...)
	allRev = append(dnaRev, rnaRev...)
	dnaIds = allIds[:len(dnaIds)]
	dnaFwd = allFwd[:len(dnaFwd)]
	dnaRev = allRev[:len(dnaRev)]
	rnaIds = allIds[len(dnaIds):]
	rnaFwd = allFwd[len(dnaFwd):]
	rnaRev = allRev[len(dnaRev):]
	return
}

func removeId(id string, allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev []string) (_allIds, _allFwd, _allRev, _dnaIds, _dnaFwd, _dnaRev, _rnaIds, _rnaFwd, _rnaRev []string) {
	var idx int
	for idx = range allIds {
		if allIds[idx] == id {
			break
		}
	}

	allIds = slices.Delete(allIds, idx, idx+1)
	allFwd = slices.Delete(allFwd, idx, idx+1)
	allRev = slices.Delete(allRev, idx, idx+1)
	dnaIds = allIds[:len(dnaIds)]
	dnaFwd = allFwd[:len(dnaFwd)]
	dnaRev = allRev[:len(dnaRev)]
	rnaIds = allIds[len(dnaIds):]
	rnaFwd = allFwd[len(dnaFwd):]
	rnaRev = allRev[len(dnaRev):]

	return allIds, allFwd, allRev, dnaIds, dnaFwd, dnaRev, rnaIds, rnaFwd, rnaRev
}
