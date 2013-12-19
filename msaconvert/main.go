package main

import (
	"flag"
	"io"
	"path"

	"github.com/TuftsBCB/io/msa"
	"github.com/TuftsBCB/seq"
	"github.com/TuftsBCB/tools/util"
)

type (
	msaReader func(io.Reader) (seq.MSA, error)
	msaWriter func(io.Writer, seq.MSA) error
	msaIO     struct {
		r msaReader
		w msaWriter
	}
)

var (
	flagInFmt  = ""
	flagOutFmt = ""

	extToFmt = map[string]string{
		"fasta": "fasta", "fa": "fasta", "fas": "fasta", "ali": "fasta",
		"sto": "stockholm",
		"a2m": "a2m",
		"a3m": "a3m",
	}
	fmtToIO = map[string]msaIO{
		"fasta":     msaIO{msa.ReadFasta, msa.WriteFasta},
		"stockholm": msaIO{msa.ReadStockholm, msa.WriteStockholm},
		"a2m":       msaIO{msa.Read, msa.WriteA2M},
		"a3m":       msaIO{msa.Read, msa.WriteA3M},
	}
)

func init() {
	flag.StringVar(&flagInFmt, "infmt", flagInFmt,
		"Force the format of the input file. Legal values are fasta, "+
			"stockholm, a2m and a3m.")
	flag.StringVar(&flagOutFmt, "outfmt", flagOutFmt,
		"Force the format of the output file. Legal values are fasta, "+
			"stockholm, a2m and a3m.")

	util.FlagParse("in-msa out-msa",
		"Convert the format of an MSA file from 'in-msa' to 'out-msa'.\n"+
			"The formats are auto detected from the file's extension, but\n"+
			"they may be forced with the 'infmt' and 'outfmt' flags.")
	util.AssertNArg(2)
}

func main() {
	in, out := util.Arg(0), util.Arg(1)
	r, w := ioFromFile(in, flagInFmt).r, ioFromFile(out, flagOutFmt).w
	inf := util.OpenFile(in)
	defer inf.Close()

	msa, err := r(inf)
	util.Assert(err, "Error parsing '%s'", in)

	outf := util.CreateFile(out)
	defer outf.Close()
	util.Assert(w(outf, msa), "Error writing '%s'", out)
}

func ioFromFile(fpath, force string) msaIO {
	var fmt string
	if len(force) > 0 {
		fmt = force
	} else {
		var ok bool
		ext := path.Ext(fpath)
		if len(ext) > 0 {
			ext = ext[1:]
		}

		fmt, ok = extToFmt[ext]
		if !ok {
			util.Fatalf("Could not detect format from extension '%s'.", ext)
		}
	}

	io, ok := fmtToIO[fmt]
	if !ok {
		util.Fatalf("BUG: Could not find converters for format '%s'.", fmt)
	}
	return io
}
