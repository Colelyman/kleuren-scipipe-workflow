package main

import (
	"flag"
	"log"
	"path/filepath"

	sp "github.com/scipipe/scipipe"
)

// Define command line arguments
// Paths to exectuables
var jellyfishPath = flag.String("jellyfish",
	"jellyfish", "the path to the jellyfish executable")
var kleurenPath = flag.String("kleuren",
	"kleuren", "the path to the kleuren exectuable")
var bftPath = flag.String("bft", "bft", "the path to the BFT exectuable")

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
	kmerPattern = "kmers.txt"
)

func parseArgs() {
	flag.Parse()
	if *kmerSize%9 != 0 {
		log.Fatal("kmerSize must be a multiple of 9")
	}
}

func getGenomePaths() []string {
	genomePaths, err := filepath.Glob("*")
	if err != nil {
		log.Fatalf("No genome files were found in the directory: %s with pattern: %s and error: %s",
			*genomeDir, *genomePattern, err)
	}
	return genomePaths
}

func main() {
	parseArgs()

	wf := sp.NewWorkflow("kleuren", 4)

	countKmers := sp.NewProc(wf, "countKmers", "")
	countKmers.SetOut("kmerCount", "*.kmers.txt")
}
