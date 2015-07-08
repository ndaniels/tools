package main

import (
	"encoding/csv"
	"encoding/gob"
	"flag"
	"runtime/pprof"

	"github.com/BurntSushi/intern"

	"github.com/TuftsBCB/io/newick"
	"github.com/ndaniels/tools/util"
)

var (
	flagThreshold = 0.097702
	flagGobIt     = ""
)

func init() {
	flag.Float64Var(&flagThreshold, "threshold", flagThreshold,
		"The threshold at which to cut the tree.")
	flag.StringVar(&flagGobIt, "gobit", flagGobIt,
		"If set, alignment distances will be cached to the file given, "+
			"then mattbench-cluster will quit.")

	util.FlagUse("cpu", "cpuprof", "verbose")
	util.FlagParse(
		"(astral-alignment-dir | alignment-distances-gob) dendrogram-tree "+
			"out-clusters.csv",
		"Where `dendrogram-tree` is a file in Newick tree format.")
	if len(flagGobIt) > 0 {
		util.AssertNArg(1)
	} else {
		util.AssertNArg(3)
	}
}

func main() {
	if len(util.FlagCpuProf) > 0 {
		f := util.CreateFile(util.FlagCpuProf)
		pprof.StartCPUProfile(f)
		defer f.Close()
		defer pprof.StopCPUProfile()
	}
	if len(flagGobIt) > 0 {
		astralDir := util.Arg(0)
		dists := readAlignmentDists(astralDir)
		enc := gob.NewEncoder(util.CreateFile(flagGobIt))
		util.Assert(enc.Encode(dists), "Could not GOB encode distances")
		return
	}

	var dists *intern.Table
	if util.IsDir(util.Arg(0)) {
		dists = readAlignmentDists(util.Arg(0))
	} else {
		dec := gob.NewDecoder(util.OpenFile(util.Arg(0)))
		util.Assert(dec.Decode(&dists), "Could not GOB decode distances")
	}

	treeFile := util.Arg(1)
	outPath := util.Arg(2)

	treeReader := newick.NewReader(util.OpenFile(treeFile))
	tree, err := treeReader.ReadTree()
	util.Assert(err, "Could not read newick tree")

	csvw := csv.NewWriter(util.CreateFile(outPath))
	clusters := treeClusters(flagThreshold, dists, tree)
	util.Assert(csvw.WriteAll(clusters))
}

// clusters corresponds to a set of lists of all labels in a subtree.
type clusters [][]string

func treeClusters(
	threshold float64,
	dists *intern.Table,
	tree *newick.Tree,
) clusters {
	if len(tree.Children) == 0 {
		if len(tree.Label) > 0 {
			return clusters{[]string{tree.Label}}
		}
		return nil
	}

	// Compare all pairs in this tree. If all are within the threshold given,
	// then add this subtree as a cluster and move on. Otherwise, dig deeper.
	within := forNode(tree, func(node1 *newick.Tree) bool {
		if len(node1.Label) == 0 {
			return true
		}
		a1 := dists.Atom(node1.Label)
		return forNode(tree, func(node2 *newick.Tree) bool {
			if len(node2.Label) == 0 {
				return true
			}
			return dists.Get(a1, dists.Atom(node2.Label)) <= threshold
		})
	})
	if within {
		cluster := make([]string, 0, 10)
		forNode(tree, func(node *newick.Tree) bool {
			if len(node.Label) == 0 {
				return true
			}
			cluster = append(cluster, node.Label)
			return true
		})
		return clusters{cluster}
	}
	clusters := make(clusters, 0, len(tree.Children))
	for i := range tree.Children {
		clusters = append(
			clusters, treeClusters(threshold, dists, &tree.Children[i])...)
	}
	return clusters
}

// forNode applies `f` to each node in pre-order. If `f` returns false, then
// all traversal stops. `forNode` returns the value of the last application
// of `f`.
func forNode(t *newick.Tree, f func(*newick.Tree) bool) bool {
	if !f(t) {
		return false
	}
	for i := range t.Children {
		if !forNode(&t.Children[i], f) {
			return false
		}
	}
	return true
}
