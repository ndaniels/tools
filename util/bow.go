package util

import (
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
	"sync"

	"github.com/TuftsBCB/fragbag"
	"github.com/TuftsBCB/fragbag/bow"
	"github.com/TuftsBCB/io/fasta"
)

// ProcessBowers is a convenient wrapper around BowerOpen that processes each
// bower value in parallel and sends the resulting BOW value on the channel
// returned. The number of goroutines spawned is equivalent to N.
//
// It is appropriate for fpaths to be the arguments given to be command line
// arguments. Each directory is recursively expanded to its files, and special
// syntax is parsed as well.
//
// If `hideProgress` is true, then a progress bar will not be emitted to
// stderr.
func ProcessBowers(
	fpaths []string,
	lib fragbag.Library,
	n int,
	hideProgress bool,
) <-chan bow.Bowed {
	if n <= 0 {
		n = 1
	}
	results := make(chan bow.Bowed, n*2)
	fpaths = AllFilesFromArgs(fpaths)

	go func() {
		var progress *Progress
		totalJobs := 0
		if !hideProgress {
			totalJobs = numJobs(fpaths)
			progress = NewProgress(totalJobs)
		}

		// We use two levels of concurrency here. The first is at the level
		// of translating files into bowers. The second is at the level of
		// computing BOWs from bowers.
		// The first level is necessary because there can a large number of
		// bower files given, where each file only produces a few BOWs.
		// The second level is necessary because there is a lot of variation
		// between the number of bowers that a single file can produce. For
		// example, while most PDB files only produce a few, some FASTA files
		// can produce millions.

		files := make(chan string, n*2)
		bs := make(chan interface{}, n*2) // channel of bowers
		wgBowers := new(sync.WaitGroup)
		wgFiles := new(sync.WaitGroup)

		// goroutines for computing BOWs from bowers
		for i := 0; i < n; i++ {
			wgBowers.Add(1)
			go func() {
				defer wgBowers.Done()
				for b := range bs {
					var bw bow.Bowed
					if fragbag.IsStructure(lib) {
						lib := lib.(fragbag.StructureLibrary)
						bw = b.(bow.StructureBower).StructureBow(lib)
					} else if fragbag.IsSequence(lib) {
						lib := lib.(fragbag.SequenceLibrary)
						bw = b.(bow.SequenceBower).SequenceBow(lib)
					} else {
						Fatalf("Unknown fragment library %T", lib)
					}
					results <- bw
				}
			}()
		}

		// goroutines for translating files into bowers
		for i := 0; i < n; i++ {
			wgFiles.Add(1)
			go func() {
				defer wgFiles.Done()

				for fpath := range files {
					var err error
					for b := range BowerOpen(fpath, lib) {
						if b.Err != nil {
							err = b.Err
						} else {
							bs <- b.Bower
						}
						if IsFasta(fpath) { // each sequence counts
							progress.JobDone(err)
						}
					}
					if IsPDB(fpath) { // PDB file only counts as one job
						progress.JobDone(err)
					}
				}
			}()
		}
		for _, fpath := range fpaths {
			files <- fpath
		}

		close(files)
		wgFiles.Wait()
		close(bs)
		wgBowers.Wait()
		progress.Close()
		close(results)
	}()
	return results
}

// BowerErr corresponds to a value that is either a Bower or an error
// indicating why a Bower value could not be constructed.
type BowerErr struct {
	Bower interface{}
	Err   error
}

// BowerOpen reads the contents of `fpath` and attempts to interpret it as a
// value (or values) that implement the `bow.Bower` interface. The list of
// `bow.Bower` values returned is guaranteed to be homogenous: they will
// either be all `bow.SequenceBower` values or `bow.StructureBower` values.
//
// The actual return value of the function is a receive-only channel of BowerErr
// values. Each BowerErr value either has the `Bower` member set or has the
// `err` field set to an error that prevented the file from being opened.
// Errors in this case are reserved for files that appear to be capable of
// producing a BOW, but were unable to be read.
//
// If the fpath given cannot be detected as a bower file, then a closed empty
// channel will be returned. A warning is also emitted to stderr.
//
// `lib` is a fragment library that is used to help interpret what kind of
// value must be in `r`. For example, if `lib` is a sequence fragment library,
// then `BowerOpen` is guaranteed to return a `Bower` value that implements the
// `bow.SequenceBower` interface.
//
// As of now, `BowerOpen` can read these types of files:
//
//	File extension                 Format    Interpretation
//	*.{ent.gz,pdb,ent}             PDB       whatever `lib` is
//	*.{fasta,fas,fasta.gz,fas.gz}  FASTA     sequence
//	everything else                error     invalid
//
// Note that special syntax for PDB file names is supported. Namely, chain
// identifiers can be appended to the end of the file name, and only that chain
// will be included in the `bow.Bower` value. Otherwise, all chains in the PDB
// entry will be returned as individual `bow.Bower` values.
//
// The format is simple and easily demonstrated by examples:
//
//	1ctf.end.gz       Chains A and B
//	1ctf.ent.gz:A     Only chain A
//	1ctf.ent.gz:B     Only chain B
//	1ctf.ent.gz:A,B   Chains A and B
//
// A secondary format is also accepted. The following are equivalent to their
// corresponding examples above:
//
//	1ctf
//	1ctfA
//	1ctfB
//	1ctf:A,B
//
// Finally, `fpath` may be the name of a PDB identifier and its file path will
// be inferred from the value of the `PDB_PATH` environment variable.
// Alternatively, `fpath` may be the name of a SCOP domain, and its
// corresponding PDB file will be inferred from the value of the
// `SCOP_PDB_PATH` environment variable.
func BowerOpen(fpath string, lib fragbag.Library) <-chan BowerErr {
	if lib == nil {
		Fatalf("Files can only be converted to Fragbag frequency vectors " +
			"if a fragment library is specified.")
	}

	bowers := make(chan BowerErr, 100)
	switch {
	case IsPDB(fpath):
		go func() {
			defer close(bowers)

			entry, chains, err := PDBOpen(fpath)
			if err != nil {
				err = fmt.Errorf("Error reading '%s': %s", fpath, err)
				bowers <- BowerErr{Err: err}
				return
			}

			if fragbag.IsStructure(lib) {
				for i := range chains {
					if !chains[i].IsProtein() {
						continue
					}

					b := bow.BowerFromChain(chains[i])
					bowers <- BowerErr{Bower: b}
				}
			} else {
				for i := range chains {
					if !chains[i].IsProtein() {
						continue
					}

					s := chains[i].AsSequence()
					if s.Len() == 0 {
						Warnf("Chain '%s:%c' does not have an amino sequence.",
							entry.IdCode, chains[i].Ident)
						continue
					}
					bowers <- BowerErr{Bower: bow.BowerFromSequence(s)}
				}
			}
		}()
		return bowers
	case IsFasta(fpath):
		go func() {
			defer close(bowers)

			r, fp, err := fastaOpen(fpath)
			if err != nil {
				err = fmt.Errorf("Error reading file: %s", err)
				bowers <- BowerErr{Err: err}
				return
			}
			defer fp.Close()

			fr := fasta.NewReader(r)
			for {
				s, err := fr.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					err = fmt.Errorf("Error reading file: %s", err)
					bowers <- BowerErr{Err: err}
					return
				}
				bowers <- BowerErr{Bower: bow.BowerFromSequence(s)}
			}
		}()
		return bowers
	}
	Warnf("I don't know how to produce a Fragbag frequency vector "+
		"from the file '%s'.", fpath)
	close(bowers)
	return bowers
}

// numJobs returns an appoximate number of Bower values from the list of files
// provided. Note that a PDB file is counted as a single value even if there
// are multiple chains in it. On the other hand, FASTA files are counted for
// each individual sequence in the file.
//
// If there is a problem reading a file, it is ignored. (Presumably the error
// will be properly dealt with later.)
func numJobs(fpaths []string) int {
	count := 0
	for _, fpath := range fpaths {
		switch {
		case IsFasta(fpath):
			func() {
				r, fp, err := fastaOpen(fpath)
				if err != nil {
					return
				}
				defer fp.Close()

				n, _ := fasta.QuickSequenceCount(r)
				count += n
			}()
		case IsPDB(fpath):
			count += 1
		default:
			count += 1 // Errors result in a single call to JobDone.
		}
	}
	return count
}

// fastaOpen tries to open a FASTA file for reading. Both an io.Reader and a
// *os.File are returned. Namely, the underlying value of the io.Reader may
// not be an *os.File (e.g., it may be a Gzip reader).
//
// The file is returned so that the caller may close it.
func fastaOpen(fpath string) (io.Reader, *os.File, error) {
	fp, err := os.Open(fpath)
	if err != nil {
		return nil, nil, err
	}

	var r io.Reader
	if strings.HasSuffix(fpath, ".gz") {
		gr, err := gzip.NewReader(fp)
		if err != nil {
			fp.Close()
			return nil, nil, err
		}
		r = gr
	} else {
		r = fp
	}
	return r, fp, nil
}
