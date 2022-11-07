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
	"log"
	"sort"
	"strconv"
)

func usage() {
	fmt.Print(
		"callSNV - Call single nucleotide variants from META-CS data processed with annotateReadFamilies.\n" +
			"Usage:\n" +
			"callSNV [options] -i input.bam -b input.bed -r reference.fasta > output.vcf\n\n")
	flag.PrintDefaults()
}

func main() {
	input := flag.String("i", "", "Input bam file. Must be indexed.")
	output := flag.String("o", "stdout", "Output VCF file.")
	bedFile := flag.String("b", "", "Input bed file with coordinates of read families, read family ID, and read counts for watson and crick strands. Generated with -bed option in annotateReadFamilies.")
	ref := flag.String("r", "", "Fasta file with reference genome used to align input bam. Must be indexed.")
	totalDepth := flag.Int("a", 4, "Minimum total depth of read family for variant consideration.")
	strandedDepth := flag.Int("s", 2, "Minimum depth of independent watson and crick strands for variant consideration")
	minMapQ := flag.Int("minMapQ", 20, "Minimum mapping quality.")
	flag.Parse()

	if *input == "" || *bedFile == "" || *ref == "" {
		usage()
		log.Fatal("ERROR: must specify bam (-i), bed (-b), and fasta (-r).")
	}

	if *strandedDepth*2 > *totalDepth {
		log.Fatal("ERROR: -s * 2 should not be larger than -a")
	}

	callSNV(*input, *output, *ref, *bedFile, uint8(*minMapQ), *totalDepth, *strandedDepth)
}

func callSNV(input, output, ref, bedFile string, minMapQ uint8, minTotalDepth, minStrandedDepth int) {
	bamReader, bamHeader := sam.OpenBam(input)
	bai := sam.ReadBai(input + ".bai")
	vcfOut := fileio.EasyCreate(output)
	faSeeker := fasta.NewSeeker(ref, "")
	bedChan := bed.GoReadToChan(bedFile)

	var err error
	var watsonDepth, crickDepth int

	for b := range bedChan {
		watsonDepth, _ = strconv.Atoi(b.Annotation[0])
		crickDepth, _ = strconv.Atoi(b.Annotation[1])
		if watsonDepth+crickDepth < minTotalDepth {
			continue
		}
		if watsonDepth < minStrandedDepth || crickDepth < minStrandedDepth {
			continue
		}
		callFamily(b, bamReader, bamHeader, faSeeker, bai, minMapQ)
	}

	err = bamReader.Close()
	exception.PanicOnErr(err)
	err = faSeeker.Close()
	exception.PanicOnErr(err)
	err = vcfOut.Close()
	exception.PanicOnErr(err)
}

func callFamily(b bed.Bed, bamReader *sam.BamReader, header sam.Header, faSeeker *fasta.Seeker, bai sam.Bai, minMapQ uint8) {
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

	watsonPiles := pileup(watsonReads, header, b.ChromStart, b.ChromEnd)
	crickPiles := pileup(crickReads, header, b.ChromStart, b.ChromEnd)

	var putativeWatsonPiles, putativeCrickPiles []sam.Pile
	var watsonPileIdx, crickPileIdx int
	var maxWatsonBase, maxCrickBase dna.Base
	var refBase []dna.Base
	for {
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
		maxWatsonBase = maxBase(watsonPiles[watsonPileIdx])
		maxCrickBase = maxBase(crickPiles[crickPileIdx])

		if maxWatsonBase != maxCrickBase {
			watsonPileIdx++
			crickPileIdx++
			continue
		}
		refBase, err = fasta.SeekByName(faSeeker, header.Chroms[watsonPiles[watsonPileIdx].RefIdx].Name, int(watsonPiles[watsonPileIdx].Pos-1), int(watsonPiles[watsonPileIdx].Pos))
		dna.AllToUpper(refBase)
		exception.PanicOnErr(err)

		if maxWatsonBase == refBase[0] {
			watsonPileIdx++
			crickPileIdx++
			continue
		}

		putativeWatsonPiles = append(putativeWatsonPiles, watsonPiles[watsonPileIdx])
		putativeCrickPiles = append(putativeCrickPiles, crickPiles[crickPileIdx])

		fmt.Println(header.Chroms[watsonPiles[watsonPileIdx].RefIdx].Name, int(watsonPiles[watsonPileIdx].Pos), dna.BasesToString(refBase), string(dna.BaseToRune(maxWatsonBase)))
		fmt.Println("WF", watsonPiles[watsonPileIdx].CountF)
		fmt.Println("WR", watsonPiles[watsonPileIdx].CountR)
		fmt.Println("CF", crickPiles[crickPileIdx].CountF)
		fmt.Println("CR", crickPiles[crickPileIdx].CountR)

		watsonPileIdx++
		crickPileIdx++
	}
	fmt.Println("Finished", b.Chrom, b.ChromStart, b.ChromEnd)
}

func pileup(reads []sam.Sam, header sam.Header, refStart, refEnd int) []sam.Pile {
	if len(reads) == 0 {
		return nil
	}

	samChan := make(chan sam.Sam, len(reads))
	for i := range reads {
		samChan <- reads[i]
	}
	close(samChan)

	ans := make([]sam.Pile, 0, refEnd-refStart)
	pileChan := sam.GoPileup(samChan, header, false, nil, nil)
	for p := range pileChan {
		ans = append(ans, p)
	}
	return ans
}

func maxBase(p sam.Pile) dna.Base {
	var max dna.Base = 0
	for i := 1; i < len(p.CountF); i++ {
		if p.CountF[i]+p.CountR[i] > p.CountF[max]+p.CountR[max] {
			max = dna.Base(i)
		}
	}
	return max
}
