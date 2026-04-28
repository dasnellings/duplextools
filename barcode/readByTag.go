package barcode

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strconv"
	"strings"
)

func SplitBamByTag(in <-chan sam.Sam, tag string) <-chan (<-chan sam.Sam) {
	out := make(chan (<-chan sam.Sam), 100)
	go goSplitBamByTag(tag, in, out)
	return out
}

func goSplitBamByTag(splitTag string, in <-chan sam.Sam, out chan<- (<-chan sam.Sam)) {
	tagMap := make(map[string]chan sam.Sam)

	var qResult interface{}
	var tagValue string
	var err error
	var found bool
	var currChan chan sam.Sam
	for r := range in {
		qResult, found, err = sam.QueryTag(r, splitTag)
		exception.PanicOnErr(err)
		if !found {
			continue
		}
		tagValue = getString(qResult)
		if currChan, found = tagMap[tagValue]; found {
			currChan <- r
			continue
		}

		// did not find existing channel for identified tag value
		currChan = make(chan sam.Sam, 10)
		tagMap[tagValue] = currChan
		currChan <- r
	}

	// close all channels
	for _, c := range tagMap {
		close(c)
	}
	close(out)
}

func getString(q interface{}) string {
	//uint8, int8, uint16, int16, uint32, int32, float32, string, or a slice of interface{}
	switch v := q.(type) {
	case uint8:
		return strconv.Itoa(int(v))
	case uint16:
		return strconv.Itoa(int(v))
	case uint32:
		return strconv.Itoa(int(v))
	case int8:
		return strconv.Itoa(int(v))
	case int16:
		return strconv.Itoa(int(v))
	case int32:
		return strconv.Itoa(int(v))
	case float32:
		return strconv.FormatFloat(float64(v), 'g', 6, 32)
	case string:
		return v
	case []interface{}:
		s := new(strings.Builder)
		for i := range v {
			if i != 0 {
				s.WriteString(",")
			}
			s.WriteString(getString(v[i]))
		}
		return s.String()
	default:
		log.Panicf("ERROR: unexpected type %v", v)
		return ""
	}
}
