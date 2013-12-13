package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"path"

	"github.com/TuftsBCB/fragbag"
	"github.com/TuftsBCB/io/pdb"
	"github.com/TuftsBCB/structure"
	"github.com/TuftsBCB/tools/util"
)

var flagOverwrite = false

func init() {
	flag.BoolVar(&flagOverwrite, "overwrite", flagOverwrite,
		"When set, any existing database will be completely overwritten.")

	util.FlagParse("kolodny-brk-file struct-lib-outfile", "")
	util.AssertLeastNArg(2)
}

func main() {
	brkFile := util.Arg(0)
	saveto := util.Arg(1)

	util.AssertOverwritable(saveto, flagOverwrite)

	fbrk := util.OpenFile(brkFile)
	defer fbrk.Close()

	brkContents, err := ioutil.ReadAll(fbrk)
	util.Assert(err, "Could not read '%s'", brkFile)

	pdbFragments := bytes.Split(brkContents, []byte("TER"))
	fragments := make([][]structure.Coords, 0)
	for i, pdbFrag := range pdbFragments {
		pdbFrag = bytes.TrimSpace(pdbFrag)
		if len(pdbFrag) == 0 {
			continue
		}
		fragments = append(fragments, coords(i, pdbFrag))
	}

	savetof := util.CreateFile(saveto)
	defer savetof.Close()

	lib, err := fragbag.NewStructureAtoms(path.Base(brkFile), fragments)
	util.Assert(err)
	fragbag.Save(savetof, lib)
}

func coords(num int, atomRecords []byte) []structure.Coords {
	r := bytes.NewReader(atomRecords)
	name := fmt.Sprintf("fragment %d", num)

	entry, err := pdb.Read(r, name)
	util.Assert(err, "Fragment contents could not be read in PDB format")

	atoms := entry.OneChain().CaAtoms()
	if len(atoms) == 0 {
		util.Fatalf("Fragment %d has no ATOM coordinates.", num)
	}
	return atoms
}
