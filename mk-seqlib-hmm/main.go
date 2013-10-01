package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	path "path/filepath"
	"runtime/pprof"
	"strings"
	"sync"

	"github.com/TuftsBCB/apps/hhsuite"
	"github.com/TuftsBCB/fragbag"
	iomsa "github.com/TuftsBCB/io/msa"
	"github.com/TuftsBCB/io/pdb"
	"github.com/TuftsBCB/io/pdb/slct"
	"github.com/TuftsBCB/seq"
	"github.com/TuftsBCB/tools/util"
)

var (
	flagOverwrite = false
	flagPdbSelect = false
)

var (
	// There are two concurrent aspects going on here:
	// 1) processing entire PDB chains
	// 2) adding each part of each chain to a sequence fragment.
	// So we use two waitgroups: one for synchronizing on finishing
	// (1) and the other for synchronizing on finishing (2).
	wgPDBChains    = new(sync.WaitGroup)
	wgSeqFragments = new(sync.WaitGroup)

	// The structure library supplied by the user.
	structLib fragbag.StructureLibrary
)

func init() {
	flag.BoolVar(&flagOverwrite, "overwrite", flagOverwrite,
		"When set, any existing database will be completely overwritten.")
	flag.BoolVar(&flagPdbSelect, "pdb-select", flagPdbSelect,
		"When set, the protein list will be read as a PDB Select file.")

	util.FlagUse("cpu", "cpuprof", "verbose")
	util.FlagParse(
		"struct-lib protein-list seq-lib-outfile",
		"Where 'protein-list' is a plain text file with PDB chain\n"+
			"identifiers on each line. e.g., '1P9GA'.")
	util.AssertLeastNArg(3)
}

func main() {
	if len(util.FlagCpuProf) > 0 {
		f := util.CreateFile(util.FlagCpuProf)
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	libPath := util.Arg(0)
	protList := util.Arg(1)
	saveto := util.Arg(2)

	structLib = util.StructureLibrary(libPath)
	util.AssertOverwritable(saveto, flagOverwrite)

	tempDir, err := ioutil.TempDir("", "mk-seqlib-hmm")
	util.Assert(err, "Could not create temporary directory.")
	defer os.RemoveAll(tempDir)

	// Initialize a frequency profile for each structural fragment.
	var msas []seq.MSA
	var msaChans []chan seq.Sequence
	for i := 0; i < structLib.Size(); i++ {
		msa := seq.NewMSA()
		msa.SetLen(structLib.FragmentSize())
		msas = append(msas, msa)
		msaChans = append(msaChans, make(chan seq.Sequence))
	}

	// Now spin up a goroutine for each fragment that is responsible for
	// adding a sequence slice to itself.
	for i := 0; i < structLib.Size(); i++ {
		addToMSA(msaChans[i], &msas[i])
	}

	chainIds, numChainIds := genChains(protList)
	progressChan := progress(numChainIds)
	for i := 0; i < util.FlagCpu; i++ {
		wgPDBChains.Add(1)
		go func() {
			for chainId := range chainIds {
				progressChan <- struct{}{}

				pdbPath := util.PDBPath(chainId)
				if !util.Exists(pdbPath) {
					util.Verbosef("PDB file '%s' from chain '%s' does "+
						"not exist.", pdbPath, chainId)
					continue
				}
				_, chain := util.PDBReadId(chainId)
				structureToSequence(chain, msaChans)
			}
			wgPDBChains.Done()
		}()
	}
	wgPDBChains.Wait()
	close(progressChan)

	// We've finishing reading all the PDB inputs. Now close the channels
	// and let the sequence fragments finish.
	for i := 0; i < structLib.Size(); i++ {
		close(msaChans[i])
	}
	wgSeqFragments.Wait()

	// Finally, add the sequence fragments to a new sequence fragment
	// library and save.
	seqLib := fragbag.NewSequenceHMM(structLib.Name())
	for i := 0; i < structLib.Size(); i++ {
		fname := path.Join(tempDir, fmt.Sprintf("%d.fasta", i))
		f := util.CreateFile(fname)
		util.Assert(iomsa.WriteFasta(f, msas[i]))

		hhm, err := hhsuite.HHMakePseudo.Run(fname)
		util.Assert(err)
		util.Assert(seqLib.Add(hhm.HMM))
	}
	util.Assert(seqLib.Save(util.CreateFile(saveto)))
}

// structureToSequence uses structural fragments to categorize a segment
// of alpha-carbon atoms, and adds the corresponding residues to a
// corresponding sequence fragment.
func structureToSequence(chain *pdb.Chain, msaChans []chan seq.Sequence) {
	sequence := chain.AsSequence()
	fragSize := structLib.FragmentSize()

	// If the chain is shorter than the fragment size, we can do nothing
	// with it.
	if sequence.Len() < fragSize {
		util.Verbosef("Sequence '%s' is too short (length: %d)",
			sequence.Name, sequence.Len())
		return
	}

	limit := sequence.Len() - fragSize
	for start := 0; start <= limit; start++ {
		end := start + fragSize
		atoms := chain.SequenceCaAtomSlice(start, end)
		if atoms == nil {
			// Nothing contiguous was found (a "disordered" residue perhaps).
			// So skip this part of the chain.
			continue
		}
		bestFrag := structLib.Best(atoms)

		sliced := sequence.Slice(start, end)
		msaChans[bestFrag] <- sliced
	}
}

func addToMSA(sequences chan seq.Sequence, msa *seq.MSA) {
	wgSeqFragments.Add(1)
	go func() {
		for s := range sequences {
			// We don't use Add or AddFasta since both are
			// O(#sequences * #frag-length), which gets to be quite slow
			// in the presence of a lot of sequences.
			msa.Entries = append(msa.Entries, s)
		}
		wgSeqFragments.Done()
	}()
}

func genChains(protList string) (chan string, int) {
	ids := make([]string, 0, 100)
	file := util.OpenFile(protList)
	if flagPdbSelect {
		records, err := slct.NewReader(file).ReadAll()
		util.Assert(err)
		for _, r := range records {
			if len(r.ChainID) != 5 {
				util.Fatalf("Not a valid chain identifier: '%s'", r.ChainID)
			}
			ids = append(ids, r.ChainID)
		}
	} else {
		for _, line := range util.ReadLines(file) {
			line = strings.TrimSpace(line)
			if len(line) == 0 {
				continue
			} else if len(line) != 5 {
				util.Fatalf("Not a valid chain identifier: '%s'\n"+
					"Perhaps you forgot to set 'pdb-select'?", line)
			}
			ids = append(ids, line)
		}
	}

	// Convert chain IDs to a channel.
	// Idea: multiple goroutines can read and parse PDB files in parallel.
	chains := make(chan string)
	go func() {
		for _, id := range ids {
			chains <- id
		}
		close(chains)
	}()
	return chains, len(ids)
}

func progress(total int) chan struct{} {
	count := 0
	c := make(chan struct{})
	pf := func(ft string, v ...interface{}) { fmt.Fprintf(os.Stderr, ft, v...) }
	go func() {
		for _ = range c {
			count++
			pf("\r%d/%d (%0.2f%% complete)", count, total,
				100.0*(float64(count)/float64(total)))
		}
		pf("\n")
	}()
	return c
}
