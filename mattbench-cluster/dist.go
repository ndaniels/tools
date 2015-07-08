package main

import (
	"encoding/csv"
	"log"
	path "path/filepath"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/BurntSushi/intern"

	"github.com/ndaniels/tools/util"
)

type pair struct {
	key  [2]string
	dist float64
}

func readAlignmentDists(dir string) *intern.Table {
	dists := intern.NewTable(11000)
	threads := util.FlagCpu
	addDists := make(chan []pair)
	alignFile := make(chan string)
	done := make(chan struct{})

	go func() {
		for fileDists := range addDists {
			for _, pair := range fileDists {
				a1, a2 := dists.Atom(pair.key[0]), dists.Atom(pair.key[1])
				dists.Set(a1, a2, pair.dist)
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
