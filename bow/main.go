package main

import (
	"fmt"

	"github.com/TuftsBCB/fragbag"
	"github.com/TuftsBCB/fragbag/bow"
	"github.com/TuftsBCB/tools/util"
)

func init() {
	util.FlagUse("cpu")
	util.FlagParse("frag-lib-dir chain pdb-file out-bow",
		"Computes and outputs a BOW file for the specified chain in the\n"+
			"given PDB file. If 'out-bow' is '--', then a human readable\n"+
			"version of the BOW will be printed to stdout instead.")
	util.AssertNArg(4)
}

func main() {
	libPath := util.Arg(0)
	chain := util.Arg(1)
	pdbEntryPath := util.Arg(2)
	bowOut := util.Arg(3)

	lib := util.Library(libPath)
	structLib, ok := lib.(*fragbag.StructureLibrary)
	if !ok {
		util.Fatalf("Expected StructureLibrary but got %T instead.", lib)
	}
	entry := util.PDBRead(pdbEntryPath)

	thechain := entry.Chain(chain[0])
	if thechain == nil || !thechain.IsProtein() {
		util.Fatalf("Could not find chain with identifier '%c'.", chain[0])
	}

	bow := bow.PDBChainStructure{thechain}.StructureBOW(structLib)
	if bowOut == "--" {
		fmt.Println(bow)
	} else {
		util.BOWWrite(util.CreateFile(bowOut), bow)
	}
}
