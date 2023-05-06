package main

import (
	"flag"
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
	"strings"
	"sync"
)

func usage() {
	fmt.Print(
		"mcsBurdenCorrection - Correct META-CS mutation burden for trinucleotide (or other) context bias and allele count.\n" +
			"Applies method used for NanoSeq (for more information see Abascal et al. 2021 PMID: 33911282 \n" +
			"Usage:\n" +
			"mcsBurdenCorrection [options] -i input.vcf -b readFamilies.bed -r reference.fasta > summary.tsv\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input VCF file with final variant calls (somatic only, after filters).")
	ref := flag.String("r", "", "Reference FASTA file. Must be indexed (.fai).")
	bedfile := flag.String("b", "", "Bed file of read families used for calling. If ends were trimmed from read families for variant calling, it must also be reflected in the input bed regions. DO NOT INCLUDE FAMILIES THAT WERE EXCLUDED FROM ANALYSIS!!!")
	pad := flag.Int("pad", 1, "Number of context bases on either side of variant (e.g. 0 == T, 1 == ATG, 2 == TATGA, ...")
	output := flag.String("o", "stdout", "Output summary file.")
	verbose := flag.Int("v", 0, "Verbose output by setting to >0.")
	flag.Parse()

	if *input == "" || *ref == "" || *bedfile == "" {
		usage()
		log.Fatalln("ERROR: must have inputs for -i, -b, and -r")
	}

	mcsBurdenCorrection(*input, *bedfile, *ref, *output, *pad, *verbose)
}

func mcsBurdenCorrection(invcf, bedfile, fastafile, output string, pad, verbose int) {
	ref := fasta.NewSeeker(fastafile, fastafile+".fai")
	defer cleanup(ref)
	out := fileio.EasyCreate(output)
	defer cleanup(out)

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
	go getVcfContext(vcfContextMap, invcf, ref, pad, stepOneWg, verbose)

	// STEP 2: Get genome-wide context frequency
	var genomeContextMap map[string]int
	// MAP STRUCTURE:
	// key is context, value is count
	genomeContextMap = initGenomeMap(pad)
	// spawn goroutine for getGenomeContext for concurrent processing
	stepTwoWg := new(sync.WaitGroup)
	stepTwoWg.Add(1)
	go getGenomeContext(genomeContextMap, fastafile, pad, stepTwoWg, verbose)

	// STEP 3: Get experimental context frequency
	var observedContextMap map[string]int
	// MAP STRUCTURE:
	// key is context, value is count
	observedContextMap = initGenomeMap(pad) // same structure used in step 2
	// spawn goroutine for getObservedContext for concurrent processing
	stepThreeWg := new(sync.WaitGroup)
	stepThreeWg.Add(1)
	go getObservedContext(observedContextMap, bedfile, ref, pad, stepThreeWg, verbose)

	stepOneWg.Wait()
	stepTwoWg.Wait()
	stepThreeWg.Wait()

	// STEP 4: Determine the frequency of each context in the genome
	genomeFreqMap := make(map[string]float64)
	getContextFreq(genomeFreqMap, genomeContextMap)

	// STEP 5: Determine the frequency of each context in the observed data
	observedFreqMap := make(map[string]float64)
	getContextFreq(observedFreqMap, observedContextMap)

	// STEP 6: Get ratio of genome-wide to observed context frequencies for each context
	contextRatio := make(map[string]float64)
	for context := range genomeFreqMap {
		contextRatio[context] = genomeFreqMap[context] / observedFreqMap[context]
	}

	// STEP 7: Adjust mutation counts by context ratio
	var mutationBurden, adjCount float64
	adjVcfContextMap := make(map[string]map[string]float64)
	for mutation, contextMap := range vcfContextMap {
		adjVcfContextMap[mutation] = make(map[string]float64)
		for context, count := range contextMap {
			adjCount = float64(count) * contextRatio[context]
			adjVcfContextMap[mutation][context] = adjCount
			mutationBurden += adjCount
		}
	}

	// STEP 8: Calculate burden per base per genome
	var adjMutationBurden float64
	adjMutationBurden = mutationBurden / float64(sumMap(observedContextMap))

	fmt.Println(adjMutationBurden) // TODO make output files
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

func getObservedContext(m map[string]int, bedfile string, ref *fasta.Seeker, pad int, wg *sync.WaitGroup, verbose int) {
	families := bed.GoReadToChan(bedfile)
	var seq []dna.Base
	var err error
	for family := range families {
		seq, err = fasta.SeekByName(ref, family.Chrom, family.ChromStart, family.ChromEnd)
		exception.PanicOnErr(err)
		genomeContext(m, seq, pad)
	}
	wg.Done()
}

func getGenomeContext(m map[string]int, fastafile string, pad int, wg *sync.WaitGroup, verbose int) {
	ref := fasta.Read(fastafile)
	for _, chr := range ref {
		genomeContext(m, chr.Seq, pad)
	}
	wg.Done()
}

func genomeContext(m map[string]int, seq []dna.Base, pad int) {
	var s string
	var b []dna.Base
	var keyFound bool
	for i := range seq {
		if i < pad || i+pad >= len(seq) || !dna.DefineBase(seq[i]) {
			continue
		}

		s = dna.BasesToString(seq[i-pad : i+pad])

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

func getVcfContext(m map[string]map[string]int, vcfFile string, ref *fasta.Seeker, pad int, wg *sync.WaitGroup, verbose int) {
	vcfChan, _ := vcf.GoReadToChan(vcfFile)
	for v := range vcfChan {
		vcfContext(v, m, ref, pad, verbose)
	}
	mergeComplements(m)
	wg.Done()
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

func getOutput(m map[string]map[string]int) string {
	var lines []string
	var k1, k2 string
	var val int
	var subMap map[string]int
	for k1, subMap = range m {
		for k2, val = range subMap {
			lines = append(lines, fmt.Sprintf("%s\t%s\t%d", k1, k2, val))
		}
	}
	slices.Sort(lines)
	return "Variant\tContext\tCount\n" + strings.Join(lines, "\n") + "\n"
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}
