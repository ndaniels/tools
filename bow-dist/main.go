package main

import (
	"fmt"
	"math"

	"github.com/ndaniels/tools/util"
)

func init() {
	util.FlagParse("bow1 bow2", "")
	util.AssertNArg(2)
}

func main() {
	b1 := util.BowRead(util.Arg(0))
	b2 := util.BowRead(util.Arg(1))
	fmt.Printf("%0.4f\n", math.Abs(b1.Bow.Cosine(b2.Bow)))
}
