package main

import (
	"encoding/csv"
	"log"
	path "path/filepath"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/TuftsBCB/tools/util"
)

type distPairs struct {
	Dists   [][]float64
	Atoms   map[string]uint16
	Current uint16
}

type pair struct {
	key  [2]string
	dist float64
}

func (dp *distPairs) addPair(p pair) {
	k1, k2 := dp.intern(p.key[0]), dp.intern(p.key[1])
	if dp.Dists == nil {
		dp.Dists = make([][]float64, 1000)
	}
	if k1 >= uint16(len(dp.Dists)) {
		newLength := uint16(2 * len(dp.Dists))
		if k1 >= newLength {
			newLength = k1 + 1
		}
		n := make([][]float64, newLength)
		copy(n, dp.Dists)
		dp.Dists = n
	}
	if dp.Dists[k1] == nil {
		dp.Dists[k1] = make([]float64, 1000)
	}
	if k2 >= uint16(len(dp.Dists[k1])) {
		newLength := uint16(2 * len(dp.Dists[k1]))
		if k2 >= newLength {
			newLength = k2 + 1
		}
		n := make([]float64, newLength)
		copy(n, dp.Dists[k1])
		dp.Dists[k1] = n
	}
	dp.Dists[k1][k2] = p.dist
}

func (dp *distPairs) intern(s string) uint16 {
	if i, ok := dp.Atoms[s]; ok {
		return i
	}
	dp.Atoms[s] = dp.Current
	dp.Current++
	if dp.Current == 0 {
		panic("string interning overflow")
	}
	return dp.Current - 1
}

func readAlignmentDists(dir string) *distPairs {
	dists := &distPairs{
		Dists:   make([][]float64, 11000),
		Atoms:   make(map[string]uint16, 11000),
		Current: 0,
	}
	threads := util.FlagCpu
	addDists := make(chan []pair)
	alignFile := make(chan string)
	done := make(chan struct{})

	go func() {
		for fileDists := range addDists {
			for _, pair := range fileDists {
				dists.addPair(pair)
			}
		}
		done <- struct{}{}
	}()

	wg := new(sync.WaitGroup)
	for i := 0; i < threads; i++ {
		wg.Add(1)
		go func() {
			for fpath := range alignFile {
				log.Printf("Reading %s (%s)", fpath, time.Now())

				f := util.OpenFile(fpath)
				defer f.Close()

				csvr := csv.NewReader(f)
				csvr.Comma = '\t'
				csvr.TrimLeadingSpace = true
				csvr.FieldsPerRecord = -1 // data is poorly formatted

				records, err := csvr.ReadAll()
				util.Assert(err, "[%s]", fpath)

				fileDists := make([]pair, 0, 100000)
				for _, record := range records {
					if len(record) != 9 {
						continue
					}
					p := recordToDist(record)
					fileDists = append(fileDists, p)
				}
				addDists <- fileDists
			}
			wg.Done()
		}()
	}

	for _, fpath := range util.RecursiveFiles(dir) {
		if strings.HasPrefix(path.Base(fpath), ".") {
			continue
		}
		alignFile <- fpath
	}
	close(alignFile)
	wg.Wait()
	close(addDists)
	<-done
	return dists
}

func (dp *distPairs) dist(p1, p2 string) float64 {
	if p1 == p2 {
		return 0
	}

	var k1, k2 uint16
	var ok bool
	if p2 < p1 {
		p1, p2 = p2, p1
	}

	k1, ok = dp.Atoms[p1]
	if !ok {
		return -1
		// util.Fatalf("Could not find distance for pair (%s, %s).", p1, p2)
	}

	k2, ok = dp.Atoms[p2]
	if !ok {
		return -1
		// util.Fatalf("Could not find distance for pair (%s, %s).", p1, p2)
	}
	if k1 >= uint16(len(dp.Dists)) || k2 >= uint16(len(dp.Dists[k1])) {
		return -1
	}
	return dp.Dists[k1][k2]
}

func recordToDist(record []string) pair {
	namePieces := strings.SplitN(record[0], ".ent_", 2)
	if len(namePieces) != 2 {
		util.Fatalf("Invalid alignment pair: '%s'.", record[0])
	}
	p1, p2 := namePieces[0], namePieces[1]
	p2 = p2[0 : len(p2)-5]

	rf := func(i int) float64 { return readFloat(record[i]) }
	corelen, rmsd := rf(1), rf(2)
	l1, l2 := rf(7), rf(8)
	coreval := (2.0 * corelen) / (l1 + l2)

	dist := -6.04979701*(rmsd-coreval*corelen*0.155+1.6018) + 1000
	dist = 1.0 / dist
	dist *= 100.0
	if p1 < p2 {
		return pair{[2]string{p1, p2}, dist}
	}
	return pair{[2]string{p2, p1}, dist}
}

func readFloat(s string) float64 {
	num, err := strconv.ParseFloat(s, 64)
	util.Assert(err, "Expected float, but got '%s'.", s)
	return num
}
