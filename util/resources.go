package util

import (
	"compress/gzip"
	"encoding/gob"
	"fmt"
	"io"
	"os"
	path "path/filepath"
	"strconv"
	"strings"

	"github.com/TuftsBCB/fragbag"
	"github.com/TuftsBCB/fragbag/bow"
	"github.com/TuftsBCB/fragbag/bowdb"
	"github.com/TuftsBCB/hhfrag"
	"github.com/TuftsBCB/io/msa"
	"github.com/TuftsBCB/io/pdb"
	"github.com/TuftsBCB/seq"
)

func Library(fpath string) fragbag.Library {
	libPath := os.Getenv("FRAGLIB_PATH")
	if path.Dir(fpath) == "." && !Exists(fpath) && len(libPath) > 0 {
		fpath = path.Join(libPath, fpath)
		if !strings.HasSuffix(fpath, ".json") {
			fpath += ".json"
		}
	}
	lib, err := fragbag.Open(OpenFile(fpath))
	Assert(err, "Could not open fragment library '%s'", fpath)
	return lib
}

func StructureLibrary(path string) fragbag.StructureLibrary {
	lib, ok := Library(path).(fragbag.StructureLibrary)
	if !ok {
		Fatalf("%s (%T) is not a structure library.", path, lib)
	}
	return lib
}

func SequenceLibrary(path string) fragbag.SequenceLibrary {
	lib, ok := Library(path).(fragbag.SequenceLibrary)
	if !ok {
		Fatalf("%s (%T) is not a sequence library.", path, lib)
	}
	return lib
}

func MSA(path string) seq.MSA {
	if strings.HasSuffix(path, "a2m") || strings.HasSuffix(path, "a3m") {
		aligned, err := msa.Read(OpenFile(path))
		Assert(err, "Could not read MSA (a2m/a3m) from '%s'", path)
		return aligned
	}
	aligned, err := msa.ReadFasta(OpenFile(path))
	Assert(err, "Could not read MSA (fasta) from '%s'", path)
	return aligned
}

func OpenBOWDB(path string) *bowdb.DB {
	db, err := bowdb.OpenDB(path)
	Assert(err, "Could not open BOW database '%s'", path)
	return db
}

func PDBRead(path string) *pdb.Entry {
	entry, err := pdb.ReadPDB(path)
	Assert(err, "Could not open PDB file '%s'", path)
	return entry
}

// PDBPath takes a PDB identifier (e.g., "1ctf" or "1ctfA") and returns
// the full path to the PDB file on the file system.
//
// The PDB_PATH environment variable must be set.
func PDBPath(pid string) string {
	if !IsPDBID(pid) && !IsChainID(pid) {
		Fatalf("PDB ids must contain 4 or 5 characters, but '%s' has %d.",
			pid, len(pid))
	}
	pdbPath := os.Getenv("PDB_PATH")
	if len(pdbPath) == 0 || !IsDir(pdbPath) {
		Fatalf("The PDB_PATH environment variable must be set to open " +
			"PDB chains by just their ID.\n" +
			"PDB_PATH should be set to the directory containing a full " +
			"copy of the PDB database.")
	}

	pdbid := strings.ToLower(pid[0:4])
	group := pdbid[1:3]
	basename := fmt.Sprintf("pdb%s.ent.gz", pdbid)
	return path.Join(pdbPath, group, basename)
}

// ScopPath takes a SCOP identifier (e.g., "d3ciua1" or "d1g09c_") and returns
// the full path to the PDB file on the file system.
//
// The SCOP_PDB_PATH environment variable must be set.
func ScopPath(pid string) string {
	if len(pid) != 7 {
		Fatalf("SCOP domain ids must contain 7 characters, but '%s' has %d.",
			pid, len(pid))
	}
	pdbPath := os.Getenv("SCOP_PDB_PATH")
	if len(pdbPath) == 0 || !IsDir(pdbPath) {
		Fatalf("The SCOP_PDB_PATH environment variable must be set to open " +
			"PDB files of SCOP domain by just their ID.\n" +
			"SCOP_PDB_PATH should be set to the directory containing a full " +
			"copy of the SCOP database as PDB formatted files.")
	}

	group := pid[2:4]
	basename := fmt.Sprintf("%s.ent", pid)
	return path.Join(pdbPath, group, basename)
}

func PDBReadId(pid string) (*pdb.Entry, *pdb.Chain) {
	e := PDBRead(PDBPath(pid))
	if IsChainID(pid) {
		chain := e.Chain(pid[4])
		if chain == nil {
			Fatalf("Could not find chain '%s' in PDB entry '%s'.", pid[4], pid)
		}
		return e, chain
	}
	return e, nil
}

func GetFmap(fpath string) *hhfrag.FragmentMap {
	var fmap *hhfrag.FragmentMap
	var err error

	switch {
	case IsFasta(fpath):
		fmap, err = HHfragConf.MapFromFasta(FlagPdbHhmDB, FlagSeqDB, fpath)
		Assert(err, "Could not generate map from '%s'", fpath)
	case IsFmap(fpath):
		fmap = FmapRead(fpath)
	default:
		Fatalf("File '%s' is not a fasta or fmap file.", fpath)
	}

	return fmap
}

func FmapRead(path string) *hhfrag.FragmentMap {
	var fmap *hhfrag.FragmentMap
	f := OpenFile(path)
	defer f.Close()

	r := gob.NewDecoder(f)
	Assert(r.Decode(&fmap), "Could not GOB decode fragment map '%s'", path)
	return fmap
}

func FmapWrite(w io.Writer, fmap *hhfrag.FragmentMap) {
	encoder := gob.NewEncoder(w)
	Assert(encoder.Encode(fmap), "Could not GOB encode fragment map")
}

func BOWRead(path string) bow.BOW {
	var bow bow.BOW
	f := OpenFile(path)
	defer f.Close()

	r := gob.NewDecoder(f)
	Assert(r.Decode(&bow), "Could not GOB decode BOW '%s'", path)
	return bow
}

func BOWWrite(w io.Writer, bow bow.BOW) {
	encoder := gob.NewEncoder(w)
	Assert(encoder.Encode(bow), "Could not GOB encode BOW")
}

func OpenFile(path string) *os.File {
	f, err := os.Open(path)
	Assert(err, "Could not open file '%s'", path)
	return f
}

func CreateFile(path string) *os.File {
	f, err := os.Create(path)
	Assert(err, "Could not create file '%s'", path)
	return f
}

func ParseInt(str string) int {
	num, err := strconv.ParseInt(str, 10, 32)
	Assert(err, "Could not parse '%s' as an integer", str)
	return int(num)
}

func IsFasta(fpath string) bool {
	suffix := func(ext string) bool {
		return strings.HasSuffix(fpath, ext)
	}
	return suffix(".fasta") || suffix(".fas") ||
		suffix(".fasta.gz") || suffix(".fas.gz")
}

func OpenFasta(fpath string) io.Reader {
	if strings.HasSuffix(fpath, ".gz") {
		r, err := gzip.NewReader(OpenFile(fpath))
		Assert(err, "Could not open '%s'", fpath)
		return r
	}
	return OpenFile(fpath)
}

func IsFmap(fpath string) bool {
	return strings.HasSuffix(fpath, ".fmap")
}

func IsPDB(fpath string) bool {
	pieces := strings.Split(path.Base(fpath), ":")
	base := pieces[0]
	if path.Dir(fpath) == "." {
		if len(base) == 4 || (len(base) == 7 && base[0] == 'd') {
			return true
		}
	}

	suffix := func(ext string) bool {
		return strings.HasSuffix(base, ext)
	}
	return suffix(".ent.gz") || suffix(".pdb") || suffix(".ent")
}

func IsChainID(s string) bool {
	return len(s) == 5
}

func IsPDBID(s string) bool {
	return len(s) == 4
}
