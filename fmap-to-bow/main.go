package main

import (
	"github.com/TuftsBCB/tools/util"
)

func init() {
	util.FlagUse("cpu")
	util.FlagParse("frag-lib-dir fmap-file out-bow", "")
	util.AssertNArg(3)
}

func main() {
	lib := util.StructureLibrary(util.Arg(0))
	fmap := util.FmapRead(util.Arg(1))
	util.BOWWrite(util.CreateFile(util.Arg(2)), fmap.StructureBOW(lib))
}
