package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/barcode"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
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
		callFamily(b, bamReader, bamHeader, bai, minMapQ)
	}

	err = bamReader.Close()
	exception.PanicOnErr(err)
	err = faSeeker.Close()
	exception.PanicOnErr(err)
	err = vcfOut.Close()
	exception.PanicOnErr(err)
}

func callFamily(b bed.Bed, bamReader *sam.BamReader, header sam.Header, bai sam.Bai, minMapQ uint8) {
	var reads []sam.Sam
	var famId string
	var strand byte
	reads = sam.SeekBamRegion(bamReader, bai, b.Chrom, uint32(b.ChromStart), uint32(b.ChromEnd))
	fmt.Println("total reads:", len(reads))
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
		fmt.Println(reads[i].Extra)
		fmt.Println(famId, b.Name, strand)
		if strand == 'W' {
			watsonReads = append(watsonReads, reads[i])
		} else if strand == 'C' {
			crickReads = append(crickReads, reads[i])
		}
	}

	watsonPiles := pileup(watsonReads, header, b.ChromStart, b.ChromEnd)
	crickPiles := pileup(crickReads, header, b.ChromStart, b.ChromEnd)

	fmt.Println(b)
	fmt.Println(len(watsonReads))
	fmt.Println(len(crickReads))
	fmt.Println(len(watsonPiles))
	fmt.Println(len(crickPiles))
}

func pileup(reads []sam.Sam, header sam.Header, refStart, refEnd int) []sam.Pile {
	if len(reads) == 0 {
		return nil
	}

	samChan := make(chan sam.Sam, len(reads))
	for i := range reads {
		samChan <- reads[i]
	}

	ans := make([]sam.Pile, 0, refEnd-refStart)
	pileChan := sam.GoPileup(samChan, header, true, nil, nil)
	for p := range pileChan {
		ans = append(ans, p)
	}
	return ans
}
