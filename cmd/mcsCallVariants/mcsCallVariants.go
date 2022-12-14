package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/barcode"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
)

func usage() {
	fmt.Print(
		"mcsCallVariants - Call single nucleotide variants from META-CS data processed with annotateReadFamilies.\n" +
			"Usage:\n" +
			"callSNV [options] -i input.bam -b input.bed -r reference.fasta > output.vcf\n\n")
	flag.PrintDefaults()
}

func main() {
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to `file`")
	memprofile := flag.String("memprofile", "", "write memory profile to `file`")
	input := flag.String("i", "", "Input bam file. Must be indexed.")
	output := flag.String("o", "stdout", "Output VCF file.")
	bedFile := flag.String("b", "", "Input bed file with coordinates of read families, read family ID, and read counts for watson and crick strands. Generated with -bed option in annotateReadFamilies.")
	ref := flag.String("r", "", "Fasta file with reference genome used to align input bam. Must be indexed.")
	totalDepth := flag.Int("a", 4, "Minimum total depth of read family for variant consideration.")
	strandedDepth := flag.Int("s", 2, "Minimum depth of independent watson and crick strands for variant consideration")
	minMapQ := flag.Int("minMapQ", 20, "Minimum mapping quality.")
	minAf := flag.Float64("minAF", 0.51, "Minimum fraction of reads with alternate allele **Within a read family and within strand** to be considered a variant.")
	debugLevel := flag.Int("verbose", 0, "Level of verbosity in log.")
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	if *input == "" || *bedFile == "" || *ref == "" {
		usage()
		log.Fatal("ERROR: must specify bam (-i), bed (-b), and fasta (-r).")
	}

	if *strandedDepth*2 > *totalDepth {
		log.Fatal("ERROR: -s * 2 should not be larger than -a")
	}

	mcsCallVariants(*input, *output, *ref, *bedFile, uint8(*minMapQ), *totalDepth, *strandedDepth, *minAf, *debugLevel)

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}

func mcsCallVariants(input, output, ref, bedFile string, minMapQ uint8, minTotalDepth, minStrandedDepth int, minAf float64, debugLevel int) {
	bamReader, bamHeader := sam.OpenBam(input)
	bai := sam.ReadBai(input + ".bai")
	vcfOut := fileio.EasyCreate(output)
	vcf.NewWriteHeader(vcfOut, makeVcfHeader(input, ref))
	faSeeker := fasta.NewSeeker(ref, "")
	bedChan := bed.GoReadToChan(bedFile)

	var err error
	var watsonDepth, crickDepth int
	var familyVariants []vcf.Vcf

	for b := range bedChan {
		watsonDepth, _ = strconv.Atoi(b.Annotation[0])
		crickDepth, _ = strconv.Atoi(b.Annotation[1])
		if watsonDepth+crickDepth < minTotalDepth {
			continue
		}
		if watsonDepth < minStrandedDepth || crickDepth < minStrandedDepth {
			continue
		}
		familyVariants = callFamily(b, bamReader, bamHeader, faSeeker, bai, minMapQ, watsonDepth, crickDepth, minAf, minTotalDepth, minStrandedDepth, debugLevel)
		vcf.WriteVcfToFileHandle(vcfOut, familyVariants)
		if debugLevel > 1 {
			log.Println("Finished Read Family:", b)
		}
	}

	err = bamReader.Close()
	exception.PanicOnErr(err)
	err = faSeeker.Close()
	exception.PanicOnErr(err)
	err = vcfOut.Close()
	exception.PanicOnErr(err)
}

func callFamily(b bed.Bed, bamReader *sam.BamReader, header sam.Header, faSeeker *fasta.Seeker, bai sam.Bai, minMapQ uint8, expectedWatsonDepth, expectedCrickDepth int, minAf float64, minTotalDepth, minStrandedDepth int, debugLevel int) []vcf.Vcf {
	var reads []sam.Sam
	var famId string
	var strand byte
	var err error
	reads = sam.SeekBamRegion(bamReader, bai, b.Chrom, uint32(b.ChromStart), uint32(b.ChromEnd))
	watsonReads := make([]sam.Sam, 0, len(reads))
	crickReads := make([]sam.Sam, 0, len(reads))
	for i := range reads {
		if reads[i].MapQ < minMapQ {
			continue
		}
		sam.ParseExtra(&reads[i])
		famId = barcode.GetRF(&reads[i])
		if famId != b.Name {
			continue
		}
		strand = barcode.GetRS(&reads[i])
		if strand == 'W' {
			watsonReads = append(watsonReads, reads[i])
		} else if strand == 'C' {
			crickReads = append(crickReads, reads[i])
		}
	}

	sort.Slice(watsonReads, func(i, j int) bool {
		return watsonReads[i].Pos < watsonReads[j].Pos
	})
	sort.Slice(crickReads, func(i, j int) bool {
		return crickReads[i].Pos < crickReads[j].Pos
	})

	watsonPiles := pileup(watsonReads, header)
	crickPiles := pileup(crickReads, header)

	if debugLevel > 0 && (len(watsonReads) != expectedWatsonDepth || len(crickReads) != expectedCrickDepth) {
		log.Printf("WARNING: mismatch in expected (%d/%d) and actual (%d/%d) number of reads, may be supplementary alignments were removed at\n%s\n", expectedWatsonDepth, expectedCrickDepth, len(watsonReads), len(crickReads), b)
	}

	// remove piles that fall outside the consensus start/end of the read families
	watsonPiles, crickPiles = removePositionalOutliers(watsonPiles, crickPiles, watsonReads, crickReads)

	var variants []vcf.Vcf
	var watsonPileIdx, crickPileIdx, watsonDelLen, crickDelLen int
	var watsonInsSeq, crickInsSeq, chr string
	var maxWatsonBase, maxCrickBase dna.Base
	var refBase []dna.Base
	var watsonVarType, crickVarType variantType
	var watsonAltAlleleCount, crickAltAlleleCount, watsonInsAlleleCount, crickInsAlleleCount int
	for { // pos matching between slices of watson and crick piles
		if watsonPileIdx == len(watsonPiles) || crickPileIdx == len(crickPiles) {
			break
		}
		if crickPiles[crickPileIdx].Pos > watsonPiles[watsonPileIdx].Pos {
			watsonPileIdx++
			continue
		}
		if crickPiles[crickPileIdx].Pos < watsonPiles[watsonPileIdx].Pos {
			crickPileIdx++
			continue
		}

		// matching ref position
		watsonVarType, maxWatsonBase, watsonInsSeq, watsonDelLen, watsonAltAlleleCount, watsonInsAlleleCount = maxBase(watsonPiles[watsonPileIdx])
		crickVarType, maxCrickBase, crickInsSeq, crickDelLen, crickAltAlleleCount, crickInsAlleleCount = maxBase(crickPiles[crickPileIdx])

		// special case to bias towards insertions since they are assigned to the position before the insertion
		if float64(watsonInsAlleleCount)/float64(len(watsonReads)/2) > minAf || float64(crickInsAlleleCount)/float64(len(crickReads)/2) > minAf {
			watsonVarType = insertion
			crickVarType = insertion
			watsonAltAlleleCount = watsonInsAlleleCount
			crickAltAlleleCount = crickInsAlleleCount
		}

		// exclude if watson and crick do not agree
		if watsonVarType != crickVarType {
			watsonPileIdx++
			crickPileIdx++
			continue
		}

		// exclude if watson or crick AF is less than threshold. divide by 2 for fwd/reverse reads.
		if float64(watsonAltAlleleCount)/float64(len(watsonReads)/2) < minAf || float64(crickAltAlleleCount)/float64(len(crickReads)/2) < minAf {
			watsonPileIdx++
			crickPileIdx++
			continue
		}

		// exclude if below minimum read depth
		if watsonAltAlleleCount < minStrandedDepth || crickAltAlleleCount < minStrandedDepth || watsonAltAlleleCount+crickAltAlleleCount < minTotalDepth {
			watsonPileIdx++
			crickPileIdx++
			continue
		}

		// variant-type specific filters and processing
		chr = header.Chroms[watsonPiles[watsonPileIdx].RefIdx].Name
		switch watsonVarType {
		case snv:
			if maxWatsonBase != maxCrickBase {
				watsonPileIdx++
				crickPileIdx++
				continue
			}

			refBase, err = fasta.SeekByName(faSeeker, chr, int(watsonPiles[watsonPileIdx].Pos-1), int(watsonPiles[watsonPileIdx].Pos))
			dna.AllToUpper(refBase)
			exception.PanicOnErr(err)

			if maxWatsonBase == refBase[0] {
				watsonPileIdx++
				crickPileIdx++
				continue
			}
			variants = append(variants, snvToVcf(watsonPiles[watsonPileIdx], crickPiles[crickPileIdx], chr, refBase[0], maxWatsonBase, b.Name))

		case insertion:
			if watsonInsSeq != crickInsSeq {
				watsonPileIdx++
				crickPileIdx++
				continue
			}
			variants = append(variants, insToVcf(watsonPiles[watsonPileIdx], crickPiles[crickPileIdx], chr, watsonInsSeq, faSeeker, b.Name))

		case deletion:
			if watsonDelLen != crickDelLen {
				watsonPileIdx++
				crickPileIdx++
				continue
			}
			variants = append(variants, delToVcf(watsonPiles[watsonPileIdx], crickPiles[crickPileIdx], chr, watsonDelLen, faSeeker, b.Name))
		}

		watsonPileIdx++
		crickPileIdx++
	}
	return variants
}

func pileup(reads []sam.Sam, header sam.Header) []sam.Pile {
	if len(reads) == 0 {
		return nil
	}

	samChan := make(chan sam.Sam, len(reads))
	for i := range reads {
		samChan <- reads[i]
	}
	close(samChan)

	ans := make([]sam.Pile, 0, 100)
	pileChan := sam.GoPileup(samChan, header, false, nil, nil)
	for p := range pileChan {
		ans = append(ans, p)
	}
	return ans
}

type variantType byte

const (
	snv variantType = iota
	insertion
	deletion
	none
)

func maxBase(p sam.Pile) (tp variantType, snvAltBase dna.Base, insSeq string, delLen int, altAlleleCount, maxInsCount int) {
	var maxSnvCount, maxDelCount int

	// check SNV
	for i := 0; i < len(p.CountF); i++ {
		if i == int(dna.Gap) { // deletions handled below
			continue
		}
		if p.CountF[i]+p.CountR[i] > maxSnvCount {
			snvAltBase = dna.Base(i)
			maxSnvCount = p.CountF[i] + p.CountR[i]
		}
	}

	// check Del Fwd
	for key := range p.DelCountF {
		if p.DelCountF[key]+p.DelCountR[key] > maxDelCount {
			delLen = key
			maxDelCount = p.DelCountF[key] + p.DelCountR[key]
		}
	}

	// check Del Rev
	for key := range p.DelCountR {
		if p.DelCountF[key]+p.DelCountR[key] > maxDelCount {
			delLen = key
			maxDelCount = p.DelCountF[key] + p.DelCountR[key]
		}
	}

	// check Ins Fwd
	for key := range p.InsCountF {
		if p.InsCountF[key]+p.InsCountR[key] > maxInsCount {
			insSeq = key
			maxInsCount = p.InsCountF[key] + p.InsCountR[key]
		}
	}

	// check Ins Rev
	for key := range p.InsCountR {
		if p.InsCountF[key]+p.InsCountR[key] > maxInsCount {
			insSeq = key
			maxInsCount = p.InsCountF[key] + p.InsCountR[key]
		}
	}

	// score and return winner
	if maxSnvCount > maxInsCount && maxSnvCount > maxDelCount {
		tp = snv
		altAlleleCount = maxSnvCount
		return
	}

	if maxInsCount > maxDelCount {
		tp = insertion
		altAlleleCount = maxInsCount
		return
	}

	if delLen > 0 {
		tp = deletion
		altAlleleCount = maxDelCount
		return
	}

	tp = none
	return
}

func snvToVcf(watsonPile, crickPile sam.Pile, chr string, refBase, altBase dna.Base, readFamily string) vcf.Vcf {
	var v vcf.Vcf
	v.Chr = chr
	v.Pos = int(watsonPile.Pos)
	v.Ref = string(dna.BaseToRune(refBase))
	v.Alt = []string{string(dna.BaseToRune(altBase))}
	v.Filter = "."
	v.Info = "."
	v.Id = "."
	v.Format = []string{"GT", "DP", "WS", "CS", "RF"}

	var totalDepth, watsonDepth, crickDepth string
	totalDepth = fmt.Sprint(calcDepth(watsonPile) + calcDepth(crickPile))
	watsonDepth = fmt.Sprint(watsonPile.CountF[altBase] + watsonPile.CountR[altBase])
	crickDepth = fmt.Sprint(crickPile.CountF[altBase] + crickPile.CountR[altBase])

	v.Samples = make([]vcf.Sample, 1)
	v.Samples[0].Alleles = []int16{1}
	v.Samples[0].FormatData = []string{"", totalDepth, watsonDepth, crickDepth, readFamily}
	return v
}

func insToVcf(watsonPile, crickPile sam.Pile, chr string, insSeq string, faSeeker *fasta.Seeker, readFamily string) vcf.Vcf {
	var v vcf.Vcf
	v.Chr = chr
	v.Pos = int(watsonPile.Pos)

	refBase, err := fasta.SeekByName(faSeeker, chr, int(watsonPile.Pos)-1, int(watsonPile.Pos))
	dna.AllToUpper(refBase)
	exception.PanicOnErr(err)

	v.Ref = string(dna.BaseToRune(refBase[0]))
	v.Alt = []string{string(dna.BaseToRune(refBase[0])) + insSeq}
	v.Filter = "."
	v.Info = "."
	v.Id = "."
	v.Format = []string{"GT", "DP", "WS", "CS", "RF"}

	var totalDepth, watsonDepth, crickDepth string
	totalDepth = fmt.Sprint(calcDepth(watsonPile) + calcDepth(crickPile))
	watsonDepth = fmt.Sprint(watsonPile.InsCountF[insSeq] + watsonPile.InsCountR[insSeq])
	crickDepth = fmt.Sprint(crickPile.InsCountF[insSeq] + crickPile.InsCountR[insSeq])

	v.Samples = make([]vcf.Sample, 1)
	v.Samples[0].Alleles = []int16{1}
	v.Samples[0].FormatData = []string{"", totalDepth, watsonDepth, crickDepth, readFamily}
	return v
}

func delToVcf(watsonPile, crickPile sam.Pile, chr string, delLen int, faSeeker *fasta.Seeker, readFamily string) vcf.Vcf {
	var v vcf.Vcf
	v.Chr = chr
	v.Pos = int(watsonPile.Pos) - 1

	refBase, err := fasta.SeekByName(faSeeker, chr, int(watsonPile.Pos-2), int(watsonPile.Pos-1)+delLen)
	dna.AllToUpper(refBase)
	exception.PanicOnErr(err)

	v.Ref = dna.BasesToString(refBase)
	v.Alt = []string{string(dna.BaseToRune(refBase[0]))}
	v.Filter = "."
	v.Info = "."
	v.Id = "."
	v.Format = []string{"GT", "DP", "WS", "CS", "RF"}

	var totalDepth, watsonDepth, crickDepth string
	totalDepth = fmt.Sprint(calcDepth(watsonPile) + calcDepth(crickPile))
	watsonDepth = fmt.Sprint(watsonPile.DelCountF[delLen] + watsonPile.DelCountR[delLen])
	crickDepth = fmt.Sprint(crickPile.DelCountF[delLen] + crickPile.DelCountR[delLen])

	v.Samples = make([]vcf.Sample, 1)
	v.Samples[0].Alleles = []int16{1}
	v.Samples[0].FormatData = []string{"", totalDepth, watsonDepth, crickDepth, readFamily}
	return v
}

func makeVcfHeader(infile string, referenceFile string) vcf.Header {
	var header vcf.Header
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, fmt.Sprintf("##reference=%s", referenceFile))
	header.Text = append(header.Text, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	header.Text = append(header.Text, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=WS,Number=1,Type=Integer,Description=\"Watson Strand Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=CS,Number=1,Type=Integer,Description=\"Crick Strand Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=RF,Number=1,Type=Integer,Description=\"Read Family Identifier\">")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", strings.TrimSuffix(infile, ".bam")))
	return header
}

func removePositionalOutliers(watsonPiles, crickPiles []sam.Pile, watsonReads, crickReads []sam.Sam) (filteredWatsonPiles, filteredCrickPiles []sam.Pile) {
	filteredWatsonPiles = make([]sam.Pile, 0, len(watsonPiles))
	filteredCrickPiles = make([]sam.Pile, 0, len(crickPiles))

	fwdStartMap := make(map[int]int)
	fwdEndMap := make(map[int]int)
	revStartMap := make(map[int]int)
	revEndMap := make(map[int]int)

	for i := range watsonReads {
		if sam.IsPosStrand(watsonReads[i]) {
			fwdStartMap[watsonReads[i].GetChromStart()]++
			fwdEndMap[watsonReads[i].GetChromEnd()]++
		} else {
			revStartMap[watsonReads[i].GetChromStart()]++
			revEndMap[watsonReads[i].GetChromEnd()]++
		}
	}
	for i := range crickReads {
		if sam.IsPosStrand(crickReads[i]) {
			fwdStartMap[crickReads[i].GetChromStart()]++
			fwdEndMap[crickReads[i].GetChromEnd()]++
		} else {
			revStartMap[crickReads[i].GetChromStart()]++
			revEndMap[crickReads[i].GetChromEnd()]++
		}
	}

	var fwdStart, fwdEnd, revStart, revEnd, maxCount int
	for key, val := range fwdStartMap {
		if val > maxCount {
			fwdStart = key
		}
	}
	for key, val := range fwdEndMap {
		if val > maxCount {
			fwdEnd = key
		}
	}
	for key, val := range revStartMap {
		if val > maxCount {
			revStart = key
		}
	}
	for key, val := range revEndMap {
		if val > maxCount {
			revEnd = key
		}
	}

	for i := range watsonPiles {
		if (int(watsonPiles[i].Pos) > fwdStart && int(watsonPiles[i].Pos) < fwdEnd) ||
			(int(watsonPiles[i].Pos) > revStart && int(watsonPiles[i].Pos) < revEnd) {
			filteredWatsonPiles = append(filteredWatsonPiles, watsonPiles[i])
		}
	}

	for i := range crickPiles {
		if (int(crickPiles[i].Pos) > fwdStart && int(crickPiles[i].Pos) < fwdEnd) ||
			(int(crickPiles[i].Pos) > revStart && int(crickPiles[i].Pos) < revEnd) {
			filteredCrickPiles = append(filteredCrickPiles, crickPiles[i])
		}
	}

	return
}

// calcDepth returns the number of reads in the input pile
func calcDepth(s sam.Pile) int {
	var depth int
	for i := range s.CountF {
		depth += s.CountF[i] + s.CountR[i]
	}
	for _, val := range s.InsCountF {
		depth += val
	}
	for _, val := range s.InsCountR {
		depth += val
	}
	// Note that DelCount does NOT count towards depth
	return depth
}
