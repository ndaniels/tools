package main

import (
	"os"

	"github.com/TuftsBCB/io/hhm"
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

	qhhm, err := hhm.Read(fhhm)
	util.Assert(err)

	util.Assert(hhm.Write(os.Stdout, qhhm.Slice(start, end)))
}
