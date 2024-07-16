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

// SubCommands contains all valid subcommands.
// New subcommands can be added to duplextools by adding a new entry to this array.
var SubCommands = []*subcommand{
	{"extract", runExtract, "extract strand barcodes from fastq"},
	{"pair", runPair, "annotate reads with computed read family"},
	{"call", runCall, "call variants from duplex read families"},
	{"filter", runFilter, "remove likely constitutional variants from callset"},
	{"burden", runBurden, "calculate adjusted genome-wide mutation burden"},
}

func usage() {
	s := new(strings.Builder)
	s.WriteString(
		"Program: duplextools (Tools for duplex sequencing data)\n" +
			"Version: " + version + "\n" +
			"Contact: Daniel Snellings <daniel.snellings@childrens.harvard.edu>\n" +
			"\nUsage:\tduplextools <command> [options]\n\n" +
			"Commands:\n")

	// add subcommand text via tabwriter so the columns align
	w := tabwriter.NewWriter(s, 0, 8, 5, '\t', tabwriter.AlignRight)
	for i := range SubCommands {
		fmt.Fprintf(w, "\t%s\t%s\n", SubCommands[i].name, SubCommands[i].blurb)
	}
	w.Flush()
	fmt.Print(s.String())
}

// commandMap builds a map of possible subcommands keyed on the name of the subcommand
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

	// check if first argument is a valid subcommand
	command := commandMap()[flag.Arg(0)]

	// if no command is found, print the usage and return
	if command == nil {
		flag.Usage()
		return
	}

	// if command successfully found, execute
	command()
}
