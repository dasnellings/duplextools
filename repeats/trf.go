package repeats

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"strconv"
	"strings"
)

type Record struct {
	Bin           int
	Chr           string
	Start         int
	End           int
	Name          string
	Period        int
	CopyNum       float64
	ConsensusSize int
	PerMatch      int
	PerIndel      int
	Score         int
	A             int
	C             int
	G             int
	T             int
	Entropy       float64
	Seq           []dna.Base
}

func parseLine(line string) Record {
	words := strings.Split(line, "\t")
	var err error
	var ans Record
	ans.Bin, err = strconv.Atoi(words[0])
	exception.PanicOnErr(err)
	ans.Chr = words[1]
	ans.Start, err = strconv.Atoi(words[2])
	exception.PanicOnErr(err)
	ans.End, err = strconv.Atoi(words[3])
	exception.PanicOnErr(err)
	ans.Name = words[4]
	ans.Period, err = strconv.Atoi(words[5])
	exception.PanicOnErr(err)
	ans.CopyNum, err = strconv.ParseFloat(words[6], 64)
	exception.PanicOnErr(err)
	ans.ConsensusSize, err = strconv.Atoi(words[7])
	exception.PanicOnErr(err)
	ans.PerMatch, err = strconv.Atoi(words[8])
	exception.PanicOnErr(err)
	ans.PerIndel, err = strconv.Atoi(words[9])
	exception.PanicOnErr(err)
	ans.Score, err = strconv.Atoi(words[10])
	exception.PanicOnErr(err)
	ans.A, err = strconv.Atoi(words[11])
	exception.PanicOnErr(err)
	ans.C, err = strconv.Atoi(words[12])
	exception.PanicOnErr(err)
	ans.G, err = strconv.Atoi(words[13])
	exception.PanicOnErr(err)
	ans.T, err = strconv.Atoi(words[14])
	exception.PanicOnErr(err)
	ans.Entropy, err = strconv.ParseFloat(words[15], 64)
	exception.PanicOnErr(err)
	ans.Seq = dna.StringToBases(words[16])
	return ans
}

func GoReadToChan(file string) <-chan Record {
	ans := make(chan Record, 1000)
	go readToChan(file, ans)
	return ans
}

func readToChan(file string, c chan<- Record) {
	input := fileio.EasyOpen(file)
	var line string
	var done bool
	for line, done = fileio.EasyNextRealLine(input); !done; line, done = fileio.EasyNextRealLine(input) {
		c <- parseLine(line)
	}
	close(c)
}

//func IsPerfectRepeat(r *Record) bool {
//	length := r.End - r.Start
//	if length%len(r.Seq) != 0 {
//		return false
//	}
//	repeatedUnits := length / len(r.Seq)
//	var perfNumBases [4]int
//	for i := range r.Seq {
//		perfNumBases[r.Seq[i]] += repeatedUnits
//	}
//
//	if perfNumBases[dna.A] != r.A || perfNumBases[dna.C] != r.C || perfNumBases[dna.G] != r.G || perfNumBases[dna.T] != r.T {
//		return false
//	}
//
//	return true
//}
