package main

import (
	"fmt"

	"github.com/TuftsBCB/tools/util"
)

func init() {
	util.FlagParse("seq-lib", "")
	util.AssertNArg(1)
}

func main() {
	libPath := util.Arg(0)

	seqLib := util.SequenceLibrary(libPath)
	fmt.Println(seqLib)
	for i := 0; i < seqLib.Size(); i++ {
		fmt.Printf("%s\n\n", seqLib.FragmentString(i))
	}
}
