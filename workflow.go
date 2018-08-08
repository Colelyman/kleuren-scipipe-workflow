package main

import (
	"flag"
	"fmt"
	"log"
	"path/filepath"
	"strconv"
	"strings"

	sp "github.com/scipipe/scipipe"
)

// Define command line arguments
// Paths to exectuables
var jellyfishPath = flag.String("jellyfish",
	"./jellyfish", "the path to the jellyfish executable")
var kleurenPath = flag.String("kleuren",
	"./kleuren", "the path to the kleuren exectuable")
var bftPath = flag.String("bft", "./bft", "the path to the BFT exectuable")

// Directories
var genomeDir = flag.String("genomeDir",
	"./data",
	"the path (relative or absolute) to the directory where the genomes are stored")

// File extension patterns
var genomePattern = flag.String("genomePattern",
	"fasta", "the pattern to match genome files (with no leading .)")

// Other important flags
var kmerSize = flag.Int("k", 9, "kmer size (multiple of 9)")

// Global constant variables
const (
	kmerPattern             = "kmers.txt"
	singleJellyfishHashSize = "100M"
	multiJellyfishHashSize  = "500M"
)

func parseArgs() {
	flag.Parse()
	if *kmerSize%9 != 0 {
		log.Fatal("kmerSize must be a multiple of 9")
	}
}

func getGenomePaths() []string {
	absGenomePath, err := filepath.Abs(filepath.Join(*genomeDir, "*."+*genomePattern))
	if err != nil {
		log.Fatalf("Error getting the absolute genome path with error: %s", err)
	}
	genomePaths, err := filepath.Glob(absGenomePath)
	if err != nil {
		log.Fatalf("No genome files were found in the directory: %s with pattern: %s and error: %s",
			*genomeDir, *genomePattern, err)
	}
	return genomePaths
}

func main() {
	parseArgs()

	genomePaths := getGenomePaths()

	wf := sp.NewWorkflow("kleuren", 4)
	fmt.Println("length of genomePaths: ", len(genomePaths))

	// count the kmers for each genome
	for _, genomePath := range genomePaths {
		fmt.Println(genomePath)
		countKmers := wf.NewProc("countKmers"+genomePath,
			fmt.Sprintf("%s count -m %d -s %s -o {o:jfDB} %s",
				*jellyfishPath, *kmerSize, singleJellyfishHashSize, genomePath))
		countKmers.SetOut("jfDB", genomePath+strconv.Itoa(*kmerSize)+".jf")
		dumpKmers := wf.NewProc("dumpKmers"+genomePath,
			fmt.Sprintf("%s dump -c {i:jfDB}_* -o {o:kmerCount}", *jellyfishPath))
		dumpKmers.SetOut("kmerCount",
			"{i:jfDB|%.fasta.jf}.kmers."+strconv.Itoa(*kmerSize)+".txt")

		dumpKmers.In("jfDB").From(countKmers.Out("jfDB"))
	}

	// create the super-set of all the kmers
	absGenomePath, _ := filepath.Abs(*genomeDir)
	countSuperKmers := wf.NewProc("countSuperKmers",
		fmt.Sprintf("%s count -m %d -s %s -o {o:jfDB} %s",
			*jellyfishPath, *kmerSize, multiJellyfishHashSize, strings.Join(genomePaths, " ")))
	countSuperKmers.SetOut("jfDB",
		filepath.Join(absGenomePath, "super.kmers."+strconv.Itoa(*kmerSize)+".jf"))
	dumpSuperKmers := wf.NewProc("dumpSuperKmers",
		fmt.Sprintf("%s dump -c {i:jfDB}_* -o {o:kmerCount}", *jellyfishPath))
	dumpSuperKmers.SetOut("kmerCount",
		filepath.Join(absGenomePath, "super.kmers."+strconv.Itoa(*kmerSize)+".txt"))

	dumpSuperKmers.In("jfDB").From(countSuperKmers.Out("jfDB"))

	// Run workflow
	wf.Run()
}
