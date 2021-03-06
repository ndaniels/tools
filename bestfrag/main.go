package main

import (
	"fmt"
	"log"

	"github.com/ndaniels/esfragbag"
	"github.com/TuftsBCB/io/pdb"
	"github.com/TuftsBCB/structure"
	"github.com/ndaniels/tools/util"
)

var lib fragbag.StructureLibrary

func init() {
	u := "fraglib pdb-file [ chain-id [ start stop ] ]"
	util.FlagParse(u, "")
	util.AssertLeastNArg(2)
}

func main() {
	lib = util.StructureLibrary(util.Arg(0))
	pdbEntry := util.PDBRead(util.Arg(1))

	if util.NArg() == 2 {
		for _, chain := range pdbEntry.Chains {
			atoms := chain.CaAtoms()
			bestFragsForRegion(chain, atoms, 0, len(atoms))
		}
	} else {
		chainId := util.Arg(2)
		chain := pdbEntry.Chain(chainId[0])
		if chain == nil || !chain.IsProtein() {
			util.Fatalf("Could not find protein chain with id '%c'.", chainId)
		}
		atoms := chain.CaAtoms()

		if util.NArg() == 3 {
			bestFragsForRegion(chain, atoms, 0, len(atoms))
		} else {
			if util.NArg() != 5 {
				log.Println("Both a start and end must be provided.")
				util.Usage()
			}

			s, e := util.Arg(3), util.Arg(4)
			sn, en := util.ParseInt(s)-1, util.ParseInt(e)
			if en-sn < lib.FragmentSize() {
				util.Fatalf("The range [%s, %s] specifies %d alpha-carbon "+
					"atoms while at least %d alpha-carbon atoms are required "+
					"for the given fragment library.",
					s, e, en-sn, lib.FragmentSize())
			}
			bestFragsForRegion(chain, atoms, sn, en)
		}
	}
}

func bestFragsForRegion(chain *pdb.Chain, atoms []structure.Coords, s, e int) {
	fsize := lib.FragmentSize()
	for i := s; i <= e-fsize; i++ {
		best := lib.BestStructureFragment(atoms[i : i+fsize])
		fmt.Println(chain.Entry.IdCode, string(chain.Ident), i+1, i+fsize, best)
	}
}
