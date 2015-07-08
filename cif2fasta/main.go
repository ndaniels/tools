package main

import (
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"os"
	"path"
	"strings"

	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/io/pdbx"
	"github.com/TuftsBCB/seq"
	"github.com/ndaniels/tools/util"
)

var (
	flagChain          = ""
	flagSeparateChains = false
	flagSplit          = ""
)

func init() {
	flag.StringVar(&flagChain, "chain", flagChain,
		"This may be set to one or more chain identifiers. Only amino acids "+
			"belonging to a chain specified will be included.")
	flag.StringVar(&flagSplit, "split", flagSplit,
		"When set, each FASTA entry produced will be written to a file in the "+
			"specified directory with the PDB id code and chain identifier as "+
			"the name.")

	util.FlagParse("in-pdb-file [out-fasta-file]", "")

	if util.NArg() != 1 && util.NArg() != 2 {
		util.Usage()
	}
}

func main() {
	var f io.Reader
	var err error

	f = util.OpenFile(flag.Arg(0))
	if strings.HasSuffix(flag.Arg(0), ".gz") {
		f, err = gzip.NewReader(f)
		util.Assert(err)
	}
	cifEntry, err := pdbx.Read(f)
	util.Assert(err, "Could not read PDBx/mmCIF file")

	fasEntries := make([]seq.Sequence, 0, 5)
	for _, ent := range cifEntry.Entities {
		for _, chain := range ent.Chains {
			if !isChainUsable(chain) || len(ent.Seq) == 0 {
				continue
			}

			fasEntry := seq.Sequence{
				Name:     chainHeader(chain),
				Residues: ent.Seq,
			}
			fasEntries = append(fasEntries, fasEntry)
		}
	}
	if len(fasEntries) == 0 {
		util.Fatalf("Could not find any chains with amino acids.")
	}

	var fasOut io.Writer
	if flag.NArg() == 1 {
		fasOut = os.Stdout
	} else {
		if len(flagSplit) > 0 {
			util.Fatalf("The '--split' option is incompatible with a single " +
				"output file.")
		}
		fasOut = util.CreateFile(util.Arg(1))
	}

	if len(flagSplit) == 0 {
		util.Assert(fasta.NewWriter(fasOut).WriteAll(fasEntries),
			"Could not write FASTA file '%s'", fasOut)
	} else {
		for _, entry := range fasEntries {
			fp := path.Join(flagSplit, fmt.Sprintf("%s.fasta", entry.Name))
			out := util.CreateFile(fp)

			w := fasta.NewWriter(out)
			util.Assert(w.Write(entry), "Could not write to '%s'", fp)
			util.Assert(w.Flush(), "Could not write to '%s'", fp)
		}
	}
}

func chainHeader(chain *pdbx.Chain) string {
	ident := chain.Id
	if ident == ' ' {
		ident = 'A'
	}
	return fmt.Sprintf("%s%c", strings.ToLower(chain.Entity.Entry.Id), ident)
}

func isChainUsable(chain *pdbx.Chain) bool {
	if len(flagChain) == 0 {
		return true
	}
	for i := 0; i < len(flagChain); i++ {
		if chain.Id == flagChain[i] {
			return true
		}
	}
	return false
}
