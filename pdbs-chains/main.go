// Command pdbs-chains reads a file in the PDB Select format and outputs the
// list of PDB identifiers in the set.
package main

import (
	"flag"
	"fmt"

	"github.com/TuftsBCB/io/pdb/slct"
	"github.com/ndaniels/tools/util"
)

var (
	flagPaths = false
)

func init() {
	flag.BoolVar(&flagPaths, "paths", flagPaths,
		"When set, the full path of each PDB chain identifier will be\n"+
			"displayed, based on the value of the PDB_PATH environment\n"+
			"variable.")

	util.FlagParse("pdb-select-file",
		"Given a file in the PDB Select format, output a list of PDB chain "+
			"identifiers (one per line).")
	util.AssertNArg(1)
}

func main() {
	pdbs := util.OpenFile(flag.Arg(0))
	defer pdbs.Close()

	entries, err := slct.NewReader(pdbs).ReadAll()
	util.Assert(err)

	for _, entry := range entries {
		if flagPaths {
			fmt.Println(util.PDBPath(entry.ChainID))
		} else {
			fmt.Println(entry.ChainID)
		}
	}
}
