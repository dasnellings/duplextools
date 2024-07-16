package main

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"strings"
)

func main() {
	f1, _ := os.Open(os.Args[1])
	f2, _ := os.Open(os.Args[2])
	selectFile := bufio.NewScanner(f1)
	targetFile := bufio.NewScanner(f2)
	outFile := fileio.EasyCreate(os.Args[3])

	var line, targetLine string
	selectFile.Scan()
	moreTargets := targetFile.Scan()
	targetLine = targetFile.Text()
	for line = selectFile.Text(); selectFile.Scan(); line = selectFile.Text() {
		for !strings.HasPrefix(targetLine, line) {
			moreTargets = targetFile.Scan()
			if !moreTargets {
				fmt.Println("premature stop")
				return
			}
			targetLine = targetFile.Text()
		}
		fmt.Fprintln(outFile, targetLine)
	}
	f1.Close()
	f2.Close()
	outFile.Close()
}
