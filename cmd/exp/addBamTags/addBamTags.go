package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"addBamTags - Add tags to a bam file.\n" +
			"Usage:\n" +
			"addBamTags [options] -i input.bam -tag AUX_TAG > output.bam\n\n")
	flag.PrintDefaults()
}

// tags is a custom type that gets filled by flag.Parse()
type tags []string

// String to satisfy flag.Value interface
func (i *tags) String() string {
	return strings.Join(*i, " ")
}

// Set to satisfy flag.Value interface
func (i *tags) Set(value string) error {
	*i = append(*i, value)
	return nil
}

func main() {
	var tagsToAdd tags
	flag.Var(&tagsToAdd, "tag", "Aux tag to add to bam file. May be declared more than once to add multiple tags.")
	input := flag.String("i", "", "Input BAM file.")
	output := flag.String("o", "stdout", "Output BAM file.")
	flag.Parse()

	if *input == "" || len(tagsToAdd) == 0 {
		usage()
		log.Fatalln("ERROR: must have inputs for -i, and -tag")
	}

	for i := range tagsToAdd {
		words := strings.Split(tagsToAdd[i], ":")
		if len(words) != 3 || len(words[0]) != 2 || len(words[1]) != 1 {
			usage()
			log.Fatalln("ERROR: tag must be in form Key:Type:Value, where Key is 2 characters and Type is 1 character")
		}
	}

	addBamTags(*input, *output, tagsToAdd)
}

func addBamTags(input, output string, tags []string) {
	inChan, header := sam.GoReadToChan(input)
	outfile := fileio.EasyCreate(output)
	defer cleanup(outfile)
	out := sam.NewBamWriter(outfile, header)
	defer cleanup(out)

	for b := range inChan {
		addTag(&b, tags)
		sam.WriteToBamFileHandle(out, b, 0)
	}
}

func addTag(s *sam.Sam, tags []string) {
	sam.ParseExtra(s)
	if s.Extra != "" {
		s.Extra += "\t"
	}
	s.Extra += strings.Join(tags, "\t")
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}
