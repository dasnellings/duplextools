package gmm

import (
	"fmt"
	"github.com/guptarohit/asciigraph"
	"math"
	"math/rand"
	"testing"
	"time"
)

func TestRunMixtureModel(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	data := []float64{16, 29, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 38, 38, 40, 40, 40, 40, 40, 40}
	//data := []float64{34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 38, 38, 40, 40, 40, 40, 40, 40}

	maxIterations := 50
	maxResets := 100
	mm := new(MixtureModel)
	var converged bool
	var iterationsRun int
	for i := 0; i < 10; i++ {
		converged, iterationsRun = RunMixtureModel(data, 2, maxIterations, maxResets, mm)
		fmt.Println(mm.Means, mm.Stdev, mm.Weights, iterationsRun, converged, mm.LogLikelihood)
		//plot(data, mm)
		//for j := range mm.Data {
		//	fmt.Println(mm.Data[j], mm.posteriors[0][j], mm.posteriors[1][j])
		//}
		//fmt.Println(mm.Weights)
		//fmt.Println(mm.LogLikelihood)
		//fmt.Println(mm.Stdev)
		//fmt.Println()
	}
}

// generate data from a normal distribution with noise
func generateData(num int, mean, std float64) []float64 {
	ans := make([]float64, num)
	for i := range ans {
		ans[i] = (rand.NormFloat64() * std) + mean
	}
	return ans
}

func gaussianHist(weight, mean, stdev float64, len int) []float64 {
	y := make([]float64, 100)
	for x := range y {
		y[x] = gaussianY(float64(x), weight, mean, stdev) * float64(len) * 2
	}
	return y
}

func gaussianY(x, a, b, c float64) float64 {
	top := math.Pow(x-b, 2)
	bot := 2 * c * c
	return a * math.Exp(-top/bot)
}

func plot(observedLengths []float64, mm *MixtureModel) {
	p := make([]float64, 100)
	for i := range observedLengths {
		p[int(observedLengths[i])]++
	}

	gaussians := make([][]float64, 2)
	gaussians[0] = gaussianHist(mm.Weights[0], mm.Means[0], mm.Stdev[0], len(observedLengths))
	gaussians[1] = gaussianHist(mm.Weights[1], mm.Means[1], mm.Stdev[1], len(observedLengths))
	fmt.Println(asciigraph.Plot(p, asciigraph.Height(10), asciigraph.Precision(0)))
	fmt.Println(asciigraph.PlotMany(gaussians, asciigraph.Precision(0), asciigraph.SeriesColors(
		asciigraph.Red,
		asciigraph.Yellow,
		asciigraph.Green,
		asciigraph.Blue,
		asciigraph.Cyan,
		asciigraph.BlueViolet,
		asciigraph.Brown,
		asciigraph.Gray,
		asciigraph.Orange,
		asciigraph.Olive,
	), asciigraph.Height(10)))
}
