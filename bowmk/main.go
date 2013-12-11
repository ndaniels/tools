package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"runtime"
	"runtime/pprof"
	"sync"

	"github.com/TuftsBCB/fragbag/bow"
	"github.com/TuftsBCB/fragbag/bowdb"
	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/io/pdb"
	"github.com/TuftsBCB/tools/util"
)

var (
	flagOverwrite = false
)

func init() {
	flag.BoolVar(&flagOverwrite, "overwrite", flagOverwrite,
		"When set, any existing database will be completely overwritten.")

	util.FlagUse("cpu", "cpuprof", "verbose")
	util.FlagParse(
		"bowdb-path frag-lib-path "+
			"(protein-dir | (protein-file [protein-file ...]))",
		"Where a protein file is a FASTA, fmap or PDB file.")
	util.AssertLeastNArg(3)
}

func main() {
	dbPath := util.Arg(0)
	libPath := util.Arg(1)
	fileArgs := flag.Args()[2:]

	if flagOverwrite {
		util.Assert(os.RemoveAll(dbPath),
			"Could not remove '%s' directory", dbPath)
	}

	db, err := bowdb.CreateDB(util.Library(libPath), dbPath)
	util.Assert(err)

	if len(util.FlagCpuProf) > 0 {
		f := util.CreateFile(util.FlagCpuProf)
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	files := util.AllFilesFromArgs(fileArgs)
	progress := util.NewProgress(numJobs(files))

	fileChan := make(chan string)
	wg := new(sync.WaitGroup)
	for i := 0; i < max(1, runtime.GOMAXPROCS(0)); i++ {
		wg.Add(1)
		go func() {
			for file := range fileChan {
				addToDB(db, file, progress)
			}
			wg.Done()
		}()
	}

	for _, file := range files {
		fileChan <- file
	}

	close(fileChan)
	wg.Wait()
	progress.Close()
	util.Assert(db.Close())
}

func addToDB(db *bowdb.DB, file string, progress *util.Progress) {
	switch {
	case util.IsFasta(file):
		fp, err := os.Open(file)
		if err != nil {
			// This is a really weird error and will cause the progress
			// meter to be off. But if we're here, it means we were able
			// to read the file at startup but not now.
			progress.JobDone(fmt.Errorf(
				"Error opening FASTA file '%s': %s", file, err))
			return
		}
		freader := fasta.NewReader(fp)
		for {
			s, err := freader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				progress.JobDone(fmt.Errorf(
					"Error reading FASTA file '%s': %s", file, err))
				continue
			}
			db.Add(bow.Sequence{s})
			progress.JobDone(nil)
		}
	case util.IsPDB(file):
		entry, err := pdb.ReadPDB(file)
		if err != nil {
			progress.JobDone(fmt.Errorf(
				"Error reading PDB file '%s': %s", file, err))
			return
		}
		for _, chain := range entry.Chains {
			if chain.IsProtein() {
				db.Add(bow.PDBChainStructure{chain})
			}
		}
		progress.JobDone(nil)
	default:
		progress.JobDone(fmt.Errorf("Unrecognized protein file: '%s'.", file))
	}
}

func numJobs(files []string) int {
	count := 0
	for _, f := range files {
		switch {
		case util.IsFasta(f):
			n, err := fasta.QuickSequenceCount(util.OpenFasta(f))
			util.Assert(err, "Could not open '%s'", f)
			count += n
		case util.IsPDB(f):
			count += 1
		default:
			count += 1 // Errors result in a single call to JobDone.
		}
	}
	return count
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
