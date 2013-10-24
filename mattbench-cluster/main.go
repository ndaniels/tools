package main

import (
	"flag"
	"fmt"

	"github.com/TuftsBCB/tools/util"
)

var (
	flagThreshold = 0.097702
)

func init() {
	flag.Float64Var(&flagThreshold, "threshold", flagThreshold,
		"The threshold at which to cut the tree.")
	util.FlagUse("cpu", "cpuprof", "verbose")
	util.FlagParse(
		"astral-alignment-dir dendrogram-tree out-clusters.csv",
		"Where `dendrogram-tree` is a file in Newick tree format.")
	util.AssertLeastNArg(3)
}

func main() {
	astralDir := util.Arg(0)
	// treeFile := util.Arg(1)
	// outPath := util.Arg(2)

	dists := readAlignmentDists(astralDir)
	fmt.Println(dists.dist("d1c04b_", "d1lxea_"))
	fmt.Println(dists.dist("d1lxea_", "d1c04b_"))
}
