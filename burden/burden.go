package burden

import (
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"sort"
	"strconv"
	"strings"
	"sync"
)

func Burden(invcf, bedfile, fastafile, output string, pad, btrim, verbose int, cacheInput, cacheOutput string) {
	out := fileio.EasyCreate(output)
	defer cleanup(out)

	var cacheOut *fileio.EasyWriter
	if cacheOutput != "" {
		cacheOut = fileio.EasyCreate(cacheOutput)
		defer cleanup(cacheOut)
	}

	// STEP 1: Get the context for each variant in the vcf file
	var vcfContextMap map[string]map[string]int
	// MAP STRUCTURE:
	// Top level key of map is mutation w/o context (e.g. "T>C")
	// Second level key is reference sequence with context (e.g. "ATG")
	// Second level value is count observed.
	// e.g. map["T>C"]["ATG"] == 1
	vcfContextMap = initVcfMap(pad)
	// spawn goroutine for getSubContext for concurrent processing
	stepOneWg := new(sync.WaitGroup)
	stepOneWg.Add(1)
	go getVcfContext(vcfContextMap, invcf, fastafile, pad, stepOneWg, verbose)
	if verbose > 0 {
		log.Println("spawned thread to get context for each variant")
	}

	// STEP 2: Get genome-wide context counts
	var genomeContextMap map[string]int
	// MAP STRUCTURE:
	// key is context, value is count
	genomeContextMap = initGenomeMap(pad)
	// spawn goroutine for getGenomeContext for concurrent processing
	stepTwoWg := new(sync.WaitGroup)
	stepTwoWg.Add(1)
	var cacheErr error
	if cacheInput != "" {
		if verbose > 0 {
			log.Println("reading cache file")
		}
		cacheErr = readCache(genomeContextMap, cacheInput)
		if cacheErr == nil && verbose > 0 {
			log.Println("read genome cache successfully")
		} else if cacheErr != nil {
			log.Println(cacheErr)
		}
		stepTwoWg.Done()
	}
	if cacheInput == "" || cacheErr != nil {
		go getGenomeContext(genomeContextMap, fastafile, pad, stepTwoWg, verbose)
		if verbose > 0 {
			log.Println("spawned thread to get genome-wide context counts")
		}
	}

	// STEP 3: Get experimental context counts
	var observedContextMap map[string]int
	// MAP STRUCTURE:
	// key is context, value is count
	observedContextMap = initGenomeMap(pad) // same structure used in step 2
	// spawn goroutine for getObservedContext for concurrent processing
	stepThreeWg := new(sync.WaitGroup)
	stepThreeWg.Add(1)
	go getObservedContext(observedContextMap, bedfile, fastafile, pad, btrim, stepThreeWg, verbose)
	if verbose > 0 {
		log.Println("spawned thread to get experimentally observed context counts")
	}

	stepOneWg.Wait()
	stepTwoWg.Wait()
	stepThreeWg.Wait()

	// STEP 4: Determine the frequency of each context in the genome
	genomeFreqMap := make(map[string]float64)
	getContextFreq(genomeFreqMap, genomeContextMap)

	// STEP 5: Determine the frequency of each context in the observed data
	observedFreqMap := make(map[string]float64)
	getContextFreq(observedFreqMap, observedContextMap)

	if verbose > 0 {
		log.Println("Genome Context Values")
		log.Println(getGenomeContextOutput(genomeContextMap, genomeFreqMap))
		log.Println()
		log.Println("Observed Context Values")
		log.Println(getObservedContextOutput(observedContextMap, observedFreqMap))
		log.Println()
	}
	if cacheOut != nil {
		fmt.Fprintln(cacheOut, getGenomeContextOutput(genomeContextMap, genomeFreqMap))
	}

	// STEP 6: Get ratio of genome-wide to observed context frequencies for each context
	contextRatio := make(map[string]float64)
	for context := range genomeFreqMap {
		contextRatio[context] = genomeFreqMap[context] / observedFreqMap[context]
		if verbose > 0 {
			log.Printf("genome:observed %s ratio = %0.3f\n", context, contextRatio[context])
		}
	}

	// STEP 7: Adjust mutation counts by context ratio
	var mutationCount int
	var mutationBurden, adjCount float64
	adjVcfContextMap := make(map[string]map[string]float64)
	for mutation, contextMap := range vcfContextMap {
		adjVcfContextMap[mutation] = make(map[string]float64)
		for context, count := range contextMap {
			mutationCount += count
			adjCount = float64(count) * contextRatio[context]
			adjVcfContextMap[mutation][context] = adjCount
			mutationBurden += adjCount
		}
	}

	// STEP 8: Calculate burden per base per genome
	var adjMutationBurden float64
	var experimentalCoverage int = sumMap(observedContextMap)
	adjMutationBurden = mutationBurden / float64(experimentalCoverage)

	if verbose > 0 {
		log.Println(getVcfContextOutput(vcfContextMap, adjVcfContextMap))
	}

	fmt.Fprintln(out, getOutput(mutationCount, mutationBurden, experimentalCoverage, adjMutationBurden, genomeContextMap, observedContextMap, genomeFreqMap, observedFreqMap, contextRatio, vcfContextMap, adjVcfContextMap))
}

func getOutput(mutationCount int, adjMutationCount float64, experimentalCoverage int, adjMutationBurden float64, genomeContextMap, observedContextMap map[string]int, genomeFreqMap, observedFreqMap, contextRatio map[string]float64, vcfContextMap map[string]map[string]int, adjVcfContextMap map[string]map[string]float64) string {
	var ans []string
	ans = append(ans, fmt.Sprintf("Mutation Count:\t%d", mutationCount))
	ans = append(ans, fmt.Sprintf("Adjusted Mutation Count:\t%0.2f", adjMutationCount))
	ans = append(ans, fmt.Sprintf("Experimental Coverage:\t%d", experimentalCoverage))
	ans = append(ans, fmt.Sprintf("Adjusted Mutation Burden:\t%0.6g", adjMutationBurden))
	ans = append(ans, "#") // newline for aesthetics
	ans = append(ans, "#") // newline for aesthetics
	ans = append(ans, "#") // newline for aesthetics
	ans = append(ans, "#ContextFrequencies#\tContext\tGenomeCount\tGenomeFreq\tObservedCount\tObservedFreq\tFreqRatio")
	contextAns := getCombContextOutput(genomeContextMap, observedContextMap, genomeFreqMap, observedFreqMap, contextRatio)
	for i := range contextAns {
		ans = append(ans, "#ContextFrequencies#\t"+contextAns[i])
	}
	ans = append(ans, "#") // newline for aesthetics
	ans = append(ans, "#") // newline for aesthetics
	ans = append(ans, "#") // newline for aesthetics
	ans = append(ans, "#MutationContexts#\tVariant\tContext\tOriginalCount\tAdjustedCount")
	variantAns := vcfContextOutput(vcfContextMap, adjVcfContextMap)
	for i := range variantAns {
		ans = append(ans, "#MutationContexts#\t"+variantAns[i])
	}

	return strings.Join(ans, "\n")
}

func getContextFreq(freq map[string]float64, count map[string]int) {
	sum := float64(sumMap(count))

	for key, val := range count {
		freq[key] = float64(val) / sum
	}
}

func sumMap(m map[string]int) int {
	var ans int
	for _, val := range m {
		ans += val
	}
	return ans
}

func getObservedContext(m map[string]int, bedfile string, fastafile string, pad, btrim int, wg *sync.WaitGroup, verbose int) {
	ref := fasta.NewSeeker(fastafile, "")
	defer cleanup(ref)

	families := bed.GoReadToChan(bedfile)
	var seq []dna.Base
	var err error
	for family := range families {
		family.ChromStart += btrim
		family.ChromEnd -= btrim
		if family.ChromStart >= family.ChromEnd {
			continue
		}
		seq, err = fasta.SeekByName(ref, family.Chrom, family.ChromStart, family.ChromEnd)
		exception.PanicOnErr(err)
		genomeContext(m, seq, pad)
	}
	wg.Done()

	if verbose > 0 {
		log.Println("finished getting experimentally observed context counts")
	}
}

func getGenomeContext(m map[string]int, fastafile string, pad int, wg *sync.WaitGroup, verbose int) {
	ref := fasta.Read(fastafile)
	for _, chr := range ref {
		genomeContext(m, chr.Seq, pad)
	}
	wg.Done()

	if verbose > 0 {
		log.Println("finished getting genome-wide context counts")
	}
}

func genomeContext(m map[string]int, seq []dna.Base, pad int) {
	var s string
	var b []dna.Base
	var keyFound bool
	for i := range seq {
		if i < pad || i+pad >= len(seq) || !dna.DefineBase(seq[i]) {
			continue
		}

		s = dna.BasesToString(seq[i-pad : i+pad+1])
		// try rev comp, else exclude if invalid bases in key
		if _, keyFound = m[s]; !keyFound {
			b = dna.StringToBases(s)
			dna.ReverseComplement(b)
			s = dna.BasesToString(b)
			if _, keyFound = m[s]; !keyFound {
				continue
			}
		}

		m[s]++
	}
}

func initGenomeMap(pad int) map[string]int {
	m := make(map[string]int)

	b := "ACGT"
	var pfs []string // pfs == possible flanking sequences
	var fullContextSeq string
	permute(b, "", pad*2, &pfs) // generate all possible flanking sequences

	// we only need contexts for pyrimidines since purines are complements
	for i := range pfs {
		// sandwich reference base between possible flanking sequence to generate all possible seqs for given core base
		fullContextSeq = pfs[i][:pad] + "C" + pfs[i][pad:]
		m[fullContextSeq] = 0

		fullContextSeq = pfs[i][:pad] + "T" + pfs[i][pad:]
		m[fullContextSeq] = 0
	}

	return m
}

func getVcfContext(m map[string]map[string]int, vcfFile string, fastafile string, pad int, wg *sync.WaitGroup, verbose int) {
	ref := fasta.NewSeeker(fastafile, "")
	defer cleanup(ref)
	vcfChan, _ := vcf.GoReadToChan(vcfFile)
	for v := range vcfChan {
		vcfContext(v, m, ref, pad, verbose)
	}
	mergeComplements(m)
	wg.Done()

	if verbose > 0 {
		log.Println("finished getting context for each variant")
	}
}

func vcfContext(v vcf.Vcf, m map[string]map[string]int, ref *fasta.Seeker, pad int, verbose int) {
	var topKey, botKey string
	var keyFound bool
	var seq []dna.Base
	var numIncluded int
	var err error

	// exclude multiallelic and indels
	if !vcf.IsBiallelic(v) || !vcf.IsSubstitution(v) || v.Pos == 1 {
		if verbose > 2 { // skip indels
			log.Printf("skipping\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	topKey = v.Ref + ">" + v.Alt[0]
	// exclude if invalid bases in key
	if _, keyFound = m[topKey]; !keyFound {
		if verbose > 1 {
			log.Printf("skipping invalid bases\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	if pad > 0 {
		seq, err = fasta.SeekByName(ref, v.Chr, (v.Pos-1)-pad, (v.Pos-1)+pad+1)
	} else {
		seq = dna.StringToBases(v.Ref)
	}
	// exclude if ref sequence does not match
	if err != nil || seq[pad] != dna.StringToBase(v.Ref) {
		if verbose > 1 {
			log.Printf("skipping error fetching seq\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	botKey = dna.BasesToString(seq)
	// exclude if invalid bases in key
	if _, keyFound = m[topKey][botKey]; !keyFound {
		if verbose > 1 {
			log.Printf("skipping invalid bases in context\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	numIncluded++
	m[topKey][botKey]++
}

// initialize map with all valid permutations of nucleotides for given pad size.
func initVcfMap(pad int) map[string]map[string]int {
	m := make(map[string]map[string]int)

	// fill in all top level keys
	b := "ACGT"
	var i, j int
	for i = 0; i < len(b); i++ {
		for j = 0; j < len(b); j++ {
			if i == j {
				continue
			}
			m[string(b[i])+">"+string(b[j])] = make(map[string]int)
		}
	}

	// fill in 2nd level keys
	var pfs []string // pfs == possible flanking sequences
	var fullContextSeq string
	permute(b, "", pad*2, &pfs) // generate all possible flanking sequences
	for key, val := range m {   // loop over keys added in previous loop
		for i := range pfs {
			fullContextSeq = pfs[i][:pad] + string(key[0]) + pfs[i][pad:] // sandwich reference base (key[0]) between possible flanking sequence to generate all possible seqs for given ref basek
			val[fullContextSeq] = 0                                       // enter context sequence in map
		}
	}
	return m
}

// permute generates all possible permutations of characters present in b of length k and stores them in ans.
func permute(b string, s string, k int, ans *[]string) {
	if k == 0 {
		*ans = append(*ans, s)
		return
	}

	for i := 0; i < len(b); i++ {
		permute(b, s+string(b[i]), k-1, ans)
	}
}

func mergeComplements(m map[string]map[string]int) {
	merge(m["C>A"], m["G>T"])
	delete(m, "G>T")
	merge(m["C>G"], m["G>C"])
	delete(m, "G>C")
	merge(m["C>T"], m["G>A"])
	delete(m, "G>A")
	merge(m["T>A"], m["A>T"])
	delete(m, "A>T")
	merge(m["T>C"], m["A>G"])
	delete(m, "A>G")
	merge(m["T>G"], m["A>C"])
	delete(m, "A>C")
}

func merge(m1, m2 map[string]int) {
	var b []dna.Base
	var keyFound bool
	for key := range m1 {
		b = dna.StringToBases(key)
		dna.ReverseComplement(b)
		if _, keyFound = m2[dna.BasesToString(b)]; !keyFound {
			log.Panicf("something went horribly wrong\n could not find key %s in map below\n%v\n", dna.BasesToString(b), m2)
		}
		m1[key] += m2[dna.BasesToString(b)]
	}
}

func vcfContextOutput(orig map[string]map[string]int, adj map[string]map[string]float64) []string {
	var lines []string
	var k1, k2 string
	var subMap map[string]int
	for k1, subMap = range orig {
		for k2 = range subMap {
			lines = append(lines, fmt.Sprintf("%s\t%s\t%d\t%.2f", k1, k2, orig[k1][k2], adj[k1][k2]))
		}
	}
	slices.Sort(lines)
	return lines
}

func getVcfContextOutput(orig map[string]map[string]int, adj map[string]map[string]float64) string {
	return "#Variant\tContext\tOriginalCount\tAdjustedCount\n" + strings.Join(vcfContextOutput(orig, adj), "\n") + "\n"
}

func getGenomeContextOutput(count map[string]int, freq map[string]float64) string {
	var lines []string
	for key, val := range count {
		lines = append(lines, fmt.Sprintf("%s\t%d\t%.6f", key, val, freq[key]))
	}
	sort.Slice(lines, func(i, j int) bool {
		switch {
		case lines[i][1] < lines[j][1]:
			return true
		case lines[i][1] > lines[j][1]:
			return false
		case lines[i][0] < lines[j][0]:
			return true
		case lines[i][0] > lines[j][0]:
			return false
		default:
			return lines[i][2] < lines[j][2]
		}
	})
	return "#Context\tCount\tFrequency\n" + strings.Join(lines, "\n") + "\n"
}

func getCombContextOutput(genomeCount, obsCount map[string]int, genomeFreq, obsFreq, contextRatio map[string]float64) []string {
	var lines []string
	for key := range genomeCount {
		lines = append(lines, fmt.Sprintf("%s\t%d\t%.6f\t%d\t%.6f\t%.3f", key, genomeCount[key], genomeFreq[key], obsCount[key], obsFreq[key], contextRatio[key]))
	}
	sort.Slice(lines, func(i, j int) bool {
		switch {
		case lines[i][1] < lines[j][1]:
			return true
		case lines[i][1] > lines[j][1]:
			return false
		case lines[i][0] < lines[j][0]:
			return true
		case lines[i][0] > lines[j][0]:
			return false
		default:
			return lines[i][2] < lines[j][2]
		}
	})
	return lines
}

func getObservedContextOutput(count map[string]int, freq map[string]float64) string {
	return getGenomeContextOutput(count, freq)
}

func readCache(m map[string]int, file string) error {
	var err error
	var done, found bool
	var line, key string
	var words []string
	var val int
	in := fileio.EasyOpen(file)
	for line, done = fileio.EasyNextRealLine(in); !done; line, done = fileio.EasyNextRealLine(in) {
		words = strings.Split(line, "\t")
		if line == "" {
			continue
		}
		if len(words) < 2 {
			return errors.New("WARNING: malformed cache file, rebuilding from genome")
		}
		key = words[0]
		val, err = strconv.Atoi(words[1])
		if err != nil {
			return errors.New(fmt.Sprintf("WARNING: failed reading cache. Could not convert %s to an integer.", words[1]))
		}
		if _, found = m[key]; !found {
			return errors.New(fmt.Sprintf("WARNING: could not read from cache. Key %s present in cache file, but not in map. Make sure the cache was made with the same -pad value as current run.", key))
		}
		m[key] = val
	}
	return nil
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}
