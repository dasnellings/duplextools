package main

import (
	"flag"
	"fmt"
	"strings"
	"text/tabwriter"
)

const version string = "0.0.1"
const gonomicsVersion string = "1.0.1"

type subcommand struct {
	name     string
	function func()
	blurb    string
}

var SubCommands = []*subcommand{
	{"extract", runExtract, "extract strand barcodes from fastq"},
	{"addreadfamilies", runAddReadFamilies, "annotate reads with computed read family"},
	{"call", runCall, "call variants from duplex read families"},
	{"filter", runFilter, "filter variants from `call`"},
	{"burden", runBurden, "calculate adjusted genome-wide mutation burden"},
}

func usage() {
	s := new(strings.Builder)
	s.WriteString(
		"Program: duplextools (Tools for duplex sequencing data)\n" +
			"Version: " + version + "\n\n" +
			"Usage:\tduplextools <command> [options]\n\n" +
			"Commands:\n")

	w := tabwriter.NewWriter(s, 0, 8, 5, '\t', tabwriter.AlignRight)
	for i := range SubCommands {
		fmt.Fprintf(w, "\t%s\t%s\n", SubCommands[i].name, SubCommands[i].blurb)
	}
	w.Flush()
	fmt.Print(s.String())
}

func commandMap() map[string]func() {
	m := make(map[string]func())
	for i := range SubCommands {
		m[SubCommands[i].name] = SubCommands[i].function
	}
	return m
}

func main() {
	flag.Usage = usage
	flag.Parse()

	command := commandMap()[flag.Arg(0)]
	if command == nil {
		flag.Usage()
		return
	}

	command()
}
