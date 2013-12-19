package main

import (
	"os"

	"github.com/TuftsBCB/io/hmm"
	"github.com/TuftsBCB/tools/util"
)

func init() {
	util.FlagParse("hhm-file start end", "")
	util.AssertNArg(3)
}

func main() {
	hhmFile := util.Arg(0)
	start := util.ParseInt(util.Arg(1))
	end := util.ParseInt(util.Arg(2))

	fhhm := util.OpenFile(hhmFile)

	qhhm, err := hmm.ReadHHM(fhhm)
	util.Assert(err)

	util.Assert(hmm.WriteHHM(os.Stdout, qhhm.Slice(start, end)))
}
