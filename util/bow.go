package util

import (
	"compress/gzip"
	"fmt"
	"io"
	"os"
	path "path/filepath"
	"strings"
	"sync"

	"github.com/TuftsBCB/fragbag"
	"github.com/TuftsBCB/fragbag/bow"
	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/io/pdb"
)

type Bowered struct {
	Id   string
	Data []byte
	Bow  bow.BOW
}

// ProcessBowers is a convenient wrapper around BowerOpen that processes each
// bower value in parallel and sends the resulting BOW value on the channel
// returned. The number of goroutines spawned is equivalent to N.
//
// If `hideProgress` is true, then a progress bar will not be emitted to
// stderr.
func ProcessBowers(
	fpaths []string,
	lib fragbag.Library,
	n int,
	hideProgress bool,
) <-chan Bowered {
	results := make(chan Bowered, 100)
	go func() {
		totalJobs := 0
		if !hideProgress {
			totalJobs = numJobs(fpaths)
		}

		bs := make(chan bow.Bower, 100)
		wg := new(sync.WaitGroup)
		if n <= 0 {
			n = 1
		}
		for i := 0; i < n; i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for b := range bs {
					var bag bow.BOW
					switch lib := lib.(type) {
					case fragbag.StructureLibrary:
						bag = b.(bow.StructureBower).StructureBOW(lib)
					case fragbag.SequenceLibrary:
						bag = b.(bow.SequenceBower).SequenceBOW(lib)
					default:
						Fatalf("Unknown fragment library %T", lib)
					}
					results <- Bowered{b.Id(), b.Data(), bag}
				}
			}()
		}

		var progress *Progress
		if !hideProgress {
			progress = NewProgress(totalJobs)
		}
		for _, fpath := range fpaths {
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
		close(bs)
		wg.Wait()
		progress.Close()
		close(results)
	}()
	return results
}

// BowerErr corresponds to a value that is either a Bower or an error
// indicating why a Bower value could not be constructed.
type BowerErr struct {
	Bower bow.Bower
	Err   error
}

// BowerOpen reads the contents of `fpath` and attempts to interpret it as a
// value (or values) that implement the `bow.Bower` interface. The list of
// `bow.Bower` values returned is guaranteed to be homogenous: they will
// either be all `bow.SequenceBower` values or `bow.StructureBower` values.
//
// The actual return value of the function is a receive-only channel BowerErr
// values. Each BowerErr value either has the `Bower` member set or has the
// `err` field set to an error that prevented the file from being opened.
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

			fp, idents := pdbNameParse(fpath)
			entry, err := pdb.ReadPDB(fp)
			if err != nil {
				err = fmt.Errorf("Error reading '%s': %s", fp, err)
				bowers <- BowerErr{Err: err}
				return
			}

			var chains []*pdb.Chain
			if len(idents) == 0 {
				chains = entry.Chains
			} else {
				chains = make([]*pdb.Chain, 0, 5)
				for _, c := range idents {
					chain := entry.Chain(c)
					if chain == nil {
						Warnf("Chain '%c' does not exist for '%s'.",
							c, entry.IdCode)
						continue
					}
					chains = append(chains, chain)
				}
			}

			if fragbag.IsStructure(lib) {
				for i := range chains {
					b := bow.PDBChainStructure{chains[i]}
					bowers <- BowerErr{Bower: b}
				}
			} else {
				for i := range chains {
					s := chains[i].AsSequence()
					if s.Len() == 0 {
						Warnf("Chain '%s:%c' does not have an amino sequence.",
							entry.IdCode, chains[i].Ident)
						continue
					}
					bowers <- BowerErr{Bower: bow.Sequence{s}}
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
				bowers <- BowerErr{Bower: bow.Sequence{s}}
			}
		}()
		return bowers
	}
	Fatalf("I don't know how to produce a Fragbag frequency vector "+
		"from the file '%s'.", fpath)
	panic("unreachable")
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
			Fatalf("I don't know how to produce a Fragbag frequency vector "+
				"from the file '%s'.", fpath)
		}
	}
	return count
}

// pdbParse returns the actual file path and a list of chain identifiers.
// When the list of chain identifiers is empty, then no specific chains were
// specified.
func pdbNameParse(fpath string) (string, []byte) {
	dir, base := path.Dir(fpath), path.Base(fpath)
	pieces := strings.Split(base, ":")
	if len(pieces) > 2 {
		Fatalf("Too many colons in PDB file path '%s'.", fpath)
	}

	var idents []byte
	base = pieces[0]
	if len(pieces) > 1 {
		chains := strings.Split(pieces[1], ",")
		idents = make([]byte, len(chains))
		for i := range chains {
			if len(chains[i]) > 1 {
				Fatalf("Chain '%s' is more than one character.", chains[i])
			}
			idents[i] = byte(chains[i][0])
		}
	}

	if dir == "." && len(base) == 4 {
		return PDBPath(base), idents
	} else if dir == "." && len(base) == 7 && base[0] == 'd' {
		return ScopPath(base), idents
	}
	return path.Join(dir, base), idents
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
