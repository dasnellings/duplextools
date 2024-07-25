package extract

import (
	"fmt"
	"github.com/dasnellings/duplexTools/barcode"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

// SAM format uses ascii offset of 33 to make everything start with individual characters
// without adding 33 you get values like spaces and newlines
const asciiOffset uint8 = 33

func Extract(r1, r2, outfile, missingBcFile string) {
	readPairs := make(chan fastq.PairedEnd, 1000)
	go fastq.PairedEndToChan(r1, r2, readPairs)

	o := fileio.EasyCreate(outfile)
	bw := sam.NewBamWriter(o, sam.GenerateHeader(nil, nil, sam.Unsorted, sam.None))

	var noBcFile *fileio.EasyWriter
	var noBcWriter *sam.BamWriter
	if missingBcFile != "" {
		noBcFile = fileio.EasyCreate(missingBcFile)
		noBcWriter = sam.NewBamWriter(noBcFile, sam.GenerateHeader(nil, nil, sam.Unsorted, sam.None))
	}

	var pair fastq.PairedEnd
	var s1, s2 sam.Sam
	s1.RName = "*"
	s2.RName = "*"
	s1.RNext = "*"
	s2.RNext = "*"
	var bcFor, bcRev, bcId, extra string
	for pair = range readPairs {
		bcFor = barcode.Extract(pair.Fwd.Seq)
		bcRev = barcode.Extract(pair.Rev.Seq)
		//TODO Fuzzy Extract & trim

		if bcFor == "*" || bcRev == "*" {
			if noBcFile != nil {
				fqToSam(&pair.Fwd, &s1, true)
				fqToSam(&pair.Rev, &s2, false)
				s1.Extra = ""
				s2.Extra = ""
				sam.WriteToBamFileHandle(noBcWriter, s1, 0)
				sam.WriteToBamFileHandle(noBcWriter, s2, 0)
			}
			continue
		}

		if bcFor > bcRev {
			bcId = bcFor + "-" + bcRev
		} else {
			bcId = bcRev + "-" + bcFor
		}

		barcode.Trim(&pair.Fwd)
		barcode.Trim(&pair.Rev)

		fqToSam(&pair.Fwd, &s1, true)
		fqToSam(&pair.Rev, &s2, false)

		extra = fmt.Sprintf("AL:Z:%s\tBC:Z:%s\tBF:Z:%s\tBR:Z:%s", bcId, bcFor+"-"+bcRev, bcFor, bcRev)
		s1.Extra = extra
		s2.Extra = extra

		sam.WriteToBamFileHandle(bw, s1, 0)
		sam.WriteToBamFileHandle(bw, s2, 0)
	}

	err := bw.Close()
	exception.PanicOnErr(err)
	err = o.Close()
	exception.PanicOnErr(err)

	if noBcFile != nil {
		err = noBcWriter.Close()
		exception.PanicOnErr(err)
		err = noBcFile.Close()
		exception.PanicOnErr(err)
	}
}

func fqToSam(fq *fastq.Fastq, sam *sam.Sam, firstInPair bool) {
	sam.QName = fq.Name
	sam.Seq = fq.Seq
	for i := range fq.Qual {
		fq.Qual[i] += asciiOffset
	}
	sam.Qual = string(fq.Qual)
	if firstInPair {
		sam.Flag = 77
	} else {
		sam.Flag = 141
	}
}
