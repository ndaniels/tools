package main

import (
	"io"
	"os"
	path "path/filepath"
	"strings"

	"github.com/TuftsBCB/io/fasta"
	"github.com/ndaniels/tools/util"
)

func init() {
	util.FlagParse("fasta-file out-dir",
		"Split a single FASTA file into a set of files for each sequence.")
	util.AssertNArg(2)
}

func main() {
	rfasta := util.OpenFasta(util.Arg(0))
	dir := util.Arg(1)
	util.Assert(os.MkdirAll(dir, 0777))

	fr := fasta.NewReader(rfasta)
	for {
		s, err := fr.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			util.Assert(err)
		}

		s.Name = strings.Fields(s.Name)[0]
		fw := util.CreateFile(path.Join(dir, s.Name+".fasta"))
		w := fasta.NewWriter(fw)
		util.Assert(w.Write(s))
		util.Assert(w.Flush())
		util.Assert(fw.Close())
	}
}
