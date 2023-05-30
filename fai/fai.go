package fai

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strconv"
	"strings"
)

// Index stores the byte offset for each fasta sequencing allowing for efficient random access.
type Index struct {
	chroms  []chrOffset    // for search by index
	nameMap map[string]int // maps chr name to index in chroms
}

// String method for Index enables easy writing with the fmt package.
func (idx Index) String() string {
	answer := new(strings.Builder)
	for i := range idx.chroms {
		answer.WriteString(idx.chroms[i].String())
		answer.WriteByte('\n')
	}
	return answer.String()
}

func (idx Index) Size(chr string) int {
	return idx.chroms[idx.nameMap[chr]].len
}

// chrOffset has offset information about each reference. Equivalent to one line of a fai file.
type chrOffset struct {
	name         string // Name of this reference sequence
	len          int    // Total length of this reference sequence, in bases
	offset       int    // Offset within the FASTA file of this sequence's first base
	basesPerLine int    // The number of bases on each line
	bytesPerLine int    // The number of bytes in each line, including the newline
}

// String method for chrOffset enables easy writing with the fmt package.
func (c chrOffset) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%d\t%d", c.name, c.len, c.offset, c.basesPerLine, c.bytesPerLine)
}

// ReadIndex reads a fai index file to an Index struct that can be used for random access.
func ReadIndex(filename string) Index {
	file := fileio.EasyOpen(filename)
	var answer Index
	var curr chrOffset
	var line string
	var col []string
	var done bool
	var err error
	for line, done = fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		col = strings.Split(line, "\t")
		if len(col) != 5 {
			log.Fatalf("ERROR: malformed index file: %s\nerror on line:\n%s\n", filename, line)
		}

		curr.name = col[0]
		curr.len, err = strconv.Atoi(col[1])
		exception.PanicOnErr(err)
		curr.offset, err = strconv.Atoi(col[2])
		exception.PanicOnErr(err)
		curr.basesPerLine, err = strconv.Atoi(col[3])
		exception.PanicOnErr(err)
		curr.bytesPerLine, err = strconv.Atoi(col[4])
		exception.PanicOnErr(err)

		answer.chroms = append(answer.chroms, curr)
	}

	err = file.Close()
	exception.PanicOnErr(err)

	answer.nameMap = make(map[string]int)
	for i := range answer.chroms {
		answer.nameMap[answer.chroms[i].name] = i
	}
	return answer
}

func IndexToVcfHeader(idx Index) string {
	ans := new(strings.Builder)
	for i := range idx.chroms {
		ans.WriteString(fmt.Sprintf("##contig=<ID=%s,length=%d>\n", idx.chroms[i].name, idx.chroms[i].len))
	}
	return ans.String()
}
