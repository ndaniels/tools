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
	if !Exists(fpath) && len(libPath) > 0 {
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
	lib := Library(path)
	libStruct, ok := lib.(fragbag.StructureLibrary)
	if !ok {
		Fatalf("%s (%T) is not a structure library.", path, lib)
	}
	return libStruct
}

func SequenceLibrary(path string) fragbag.SequenceLibrary {
	lib := Library(path)
	libSeq, ok := lib.(fragbag.SequenceLibrary)
	if !ok {
		Fatalf("%s (%T) is not a sequence library.", path, lib)
	}
	return libSeq
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

func OpenBowDB(path string) *bowdb.DB {
	db, err := bowdb.Open(path)
	Assert(err, "Could not open BOW database '%s'", path)
	return db
}

func PDBOpenMust(fpath string) (*pdb.Entry, []*pdb.Chain) {
	entry, chains, err := PDBOpen(fpath)
	Assert(err)
	return entry, chains
}

func PDBOpen(fpath string) (*pdb.Entry, []*pdb.Chain, error) {
	pdbNameParse := func(fpath string) (string, []byte, string) {
		dir, base := path.Dir(fpath), path.Base(fpath)
		pieces := strings.Split(base, ":")

		var idents []byte
		base = pieces[0]
		if len(pieces) > 2 {
			Fatalf("Too many colons in PDB file path '%s'.", fpath)
		} else if len(pieces) == 2 {
			chains := strings.Split(pieces[1], ",")
			idents = make([]byte, len(chains))
			for i := range chains {
				if len(chains[i]) > 1 {
					Fatalf("Chain '%s' is more than one character.", chains[i])
				}
				idents[i] = byte(chains[i][0])
			}
		} else if len(base) == 5 { // special case for '{pdb-id}{chain-id}'
			idents = []byte{base[4]}
			base = base[0:4]
		}

		if dir == "." {
			switch len(base) {
			case 4:
				return PDBPath(base), idents, base
			case 6:
				return CathPath(base), idents, base
			case 7:
				if base[0] == 'd' {
					return ScopPath(base), idents, base
				} else {
					return CathPath(base), idents, base
				}
			}
		}
		return path.Join(dir, base), idents, ""
	}

	fp, idents, idcode := pdbNameParse(fpath)
	entry, err := pdb.ReadPDB(fp)
	if err != nil {
		err = fmt.Errorf("Error reading '%s': %s", fp, err)
		return nil, nil, err
	}
	if len(idcode) > 0 {
		if len(idcode) == 6 || (len(idcode) == 7 && idcode[0] != 'd') {
			entry.Cath = idcode
		} else if len(idcode) == 7 && idcode[0] == 'd' {
			entry.Scop = idcode
		}
	}

	var chains []*pdb.Chain
	if len(idents) == 0 {
		chains = entry.Chains
	} else {
		chains = make([]*pdb.Chain, 0, 5)
		for _, c := range idents {
			chain := entry.Chain(c)
			if chain == nil {
				Warnf("Chain '%c' does not exist for '%s'.",
					c, entry.IdCode)
				continue
			}
			chains = append(chains, chain)
		}
	}
	return entry, chains, nil
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

// CathPath takes a CATH identifier (e.g., "2h5xB03") and returns
// the full path to the PDB file on the file system.
//
// The CATH_PDB_PATH environment variable must be set.
func CathPath(pid string) string {
	if len(pid) < 6 || len(pid) > 7 {
		Fatalf("CATH domain ids must contain 6 or 7 characters, but '%s' "+
			"has %d.", pid, len(pid))
	}
	pdbPath := os.Getenv("CATH_PDB_PATH")
	if len(pdbPath) == 0 || !IsDir(pdbPath) {
		Fatalf("The CATH_PDB_PATH environment variable must be set to open " +
			"PDB files of CATH domain by just their ID.\n" +
			"CATH_PDB_PATH should be set to the directory containing a full " +
			"copy of the CATH PDB database as PDB formatted files.")
	}

	// We have to deal with some old data sets using 6-character domain IDs.
	// This is a nightmare because there doesn't appear to be an easy
	// deterministic mapping.
	if len(pid) == 6 {
		if pid[4] == '0' {
			pid_ := fmt.Sprintf("%sA%s", pid[0:4], pid[4:6])
			if p := path.Join(pdbPath, pid_); Exists(p) {
				return p
			}
		}
		pid = fmt.Sprintf("%s0%c", pid[0:5], pid[5])
	}
	return path.Join(pdbPath, pid)
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

func BowRead(path string) bow.Bowed {
	var b bow.Bowed
	f := OpenFile(path)
	defer f.Close()

	r := gob.NewDecoder(f)
	Assert(r.Decode(&b), "Could not GOB decode BOW '%s'", path)
	return b
}

func BowWrite(w io.Writer, b bow.Bowed) {
	encoder := gob.NewEncoder(w)
	Assert(encoder.Encode(b), "Could not GOB encode BOW")
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
	if path.Dir(fpath) == "." && len(base) >= 4 && len(base) <= 7 {
		return true
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
