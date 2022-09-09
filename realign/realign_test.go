package realign

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"os/exec"
	"testing"
)

func TestRealignIndels(t *testing.T) {
	in := "testdata/bwa_input.bam"
	ref := "/Users/danielsnellings/resources/hg38.fa"
	seeker := fasta.NewSeeker(ref, "")
	reads, header := sam.GoReadToChan(in)
	out := fileio.EasyCreate("testdata/out.bam")
	bw := sam.NewBamWriter(out, header)

	output := GoRealignIndels(reads, seeker)

	for r := range output {
		sam.WriteToBamFileHandle(bw, r, 0)
	}

	seeker.Close()
	bw.Close()
	out.Close()

	cmd := exec.Command("samtools", "index", "testdata/out.bam")
	cmd.Run()
}
