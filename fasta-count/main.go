package main

import (
	"fmt"

	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/tools/util"
)

func init() {
	util.FlagParse("fasta-file",
		"Quickly count the number of sequences in a fasta file.")
	util.AssertNArg(1)
}

func main() {
	rfasta := util.OpenFasta(util.Arg(0))
	count, err := fasta.QuickSequenceCount(rfasta)
	util.Assert(err)
	fmt.Println(count)
}
