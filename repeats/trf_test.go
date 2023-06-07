package repeats

import (
	"fmt"
	"testing"
)

func TestGoReadToChan(t *testing.T) {
	testfile := "/Users/danielsnellings/Desktop/simpleRepeat.txt.gz"
	//testfile := "testdata/test.txt"
	c := GoReadToChan(testfile)
	for r := range c {
		fmt.Println(r)
	}
}
