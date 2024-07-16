package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/palette"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/text"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
	"image/color"
	"math"
)

type Plottable struct {
	Data
	grid [][]float64
	minY int
	maxY int
}

func (p Plottable) Dims() (c, r int) {
	return len(p.grid), len(p.grid[0])
}

func (p Plottable) Z(c, r int) float64 {
	return p.grid[c][r]
}

func (p Plottable) X(c int) float64 {
	return float64(c)
}

func (p Plottable) Y(r int) float64 {
	return float64(r)
}

func (p Plottable) Min() float64 {
	return -1
}

func (p Plottable) Max() float64 {
	return 1
}

func PlotData(d Data) {
	var p Plottable
	var i, j, k int
	p.minY, p.maxY = minMaxReadLen(d.splitReads)
	p.Data = d
	p.grid = make([][]float64, len(d.names))
	for i = range p.grid {
		p.grid[i] = make([]float64, p.maxY+10)
	}

	// sum up total reads for each sample
	totalReads := make([]int, len(d.names))
	for i = range d.splitReads {
		for j = range d.splitReads[i] {
			totalReads[i] += len(d.splitReads[i][j])
		}
	}

	// add read lengths to grid
	for i = range d.splitReads {
		for j = range d.splitReads[i] {
			for k = range d.splitReads[i][j] {
				if j == 0 {
					p.grid[i][d.splitReads[i][j][k]] += 1
				}
				if j == 1 {
					p.grid[i][d.splitReads[i][j][k]] -= 1
				}
			}
		}
	}

	// normalize grid to number of reads in each sample
	for i = range p.grid {
		for j = range p.grid[i] {
			if totalReads[i] > 0 {
				p.grid[i][j] /= float64(totalReads[i])
			}
		}
	}

	pal := getPalette()
	hm := plotter.NewHeatMap(p, pal)
	pl := plot.New()

	pl.Add(hm)

	pl.X.Label.Text = "Sample"
	pl.Y.Label.Text = "Repeat Length"
	pl.X.Tick.Marker = sampleTicks(d.names)
	pl.X.Tick.Label.Rotation = math.Pi / 2
	pl.X.Tick.Label.YAlign = -0.35
	pl.X.Tick.Label.XAlign = text.XRight
	pl.X.Tick.Label.Font.Size = 8
	pl.X.Tick.LineStyle = draw.LineStyle{
		Color:    color.Black,
		Width:    vg.Points(0.5),
		Dashes:   []vg.Length{vg.Millimeter * 1.4},
		DashOffs: 1.4 * vg.Millimeter,
	}

	pl.Title.Text = fmt.Sprintf("%s:%d\t%s", d.v.Chr, d.v.Pos, d.v.Id)
	pl.Title.TextStyle.Font.Size = 15
	pl.Title.Padding = 0

	err := pl.Save(30*vg.Centimeter, 30*vg.Centimeter, "/Users/danielsnellings/Desktop/heatMap.pdf")
	exception.PanicOnErr(err)
	fmt.Println("created file")
}

func minMaxReadLen(r [][][]int) (min int, max int) {
	min = 100000
	max = -1
	var i, j, curr int
	for i = range r {
		for j = range r[i] {
			for _, curr = range r[i][j] {
				if curr > max {
					max = curr
				}
				if curr < min {
					min = curr
				}
			}
		}
	}
	return
}

type sampleTicks []string

func (s sampleTicks) Ticks(min, max float64) []plot.Tick {
	var ans []plot.Tick
	for i := range s {
		if float64(i) >= min && float64(i) <= max {
			ans = append(ans, plot.Tick{Value: float64(i), Label: s[i]})
		}
	}
	return ans
}

type colors []color.Color

func (c colors) Colors() []color.Color {
	return c
}
func getPalette() palette.Palette {
	var ans colors
	for i := 255; i >= 0; i-- {
		ans = append(ans, color.RGBA{255 - uint8(i), 255 - uint8(i), 255, 255})
	}
	for i := 0; i < 256; i++ {
		ans = append(ans, color.RGBA{255, 255 - uint8(i), 255 - uint8(i), 255})
	}
	return ans
}
