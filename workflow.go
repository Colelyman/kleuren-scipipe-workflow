package main

import (
	"flag"
	"log"
	"path/filepath"
	"strconv"
	"strings"
	"time"

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
var kmerSizeString = strconv.Itoa(*kmerSize)
var numMinColors = flag.Int("n", 0, "the number of minimum colors (for kleuren bubble finding)")
var maxDepth = flag.Int("d", 30, "the maximum depth for kleuren to search for bubbles")

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

	// count the kmers for each genome
	for _, genomePath := range genomePaths {
		countKmers := wf.NewProc("countKmers"+genomePath,
			"{p:jfPath} count -m {p:k} -s {p:hash} -O -o {o:jfDB} {p:genome}")
		countKmers.InParam("jfPath").FromStr(*jellyfishPath)
		countKmers.InParam("k").FromInt(*kmerSize)
		countKmers.InParam("hash").FromStr(singleJellyfishHashSize)
		countKmers.InParam("genome").FromStr(genomePath)
		countKmers.SetOut("jfDB", "{p:genome|%.fasta}.{p:k}.jf")

		dumpKmers := wf.NewProc("dumpKmers"+genomePath,
			"{p:jfPath} dump -c {i:jfDB} -o {o:kmerCount}")
		dumpKmers.InParam("jfPath").FromStr(*jellyfishPath)
		dumpKmers.SetOut("kmerCount", "{i:jfDB|%.jf}.kmers.txt")

		dumpKmers.In("jfDB").From(countKmers.Out("jfDB"))
	}

	absGenomePath, _ := filepath.Abs(*genomeDir)

	// create the color file
	colorFile := wf.NewProc("colorFile", "ls {p:genomeDir}/*.{p:k}.kmers.txt > {o:colorFile}")
	colorFile.InParam("k").FromInt(*kmerSize)
	colorFile.InParam("genomeDir").FromStr(absGenomePath)
	colorFile.SetOut("colorFile", "{p:genomeDir}/colors.txt")

	// create the super-set of all the kmers
	countSuperKmers := wf.NewProc("countSuperKmers",
		"{p:jfPath} count -m {p:k} -s {p:hash} -O -o {o:jfDB} {p:genomes}")
	countSuperKmers.InParam("jfPath").FromStr(*jellyfishPath)
	countSuperKmers.InParam("k").FromInt(*kmerSize)
	countSuperKmers.InParam("hash").FromStr(multiJellyfishHashSize)
	countSuperKmers.InParam("genomes").FromStr(strings.Join(genomePaths, " "))
	countSuperKmers.SetOut("jfDB", filepath.Join(absGenomePath, "super.kmers.{p:k}.jf"))

	dumpSuperKmers := wf.NewProc("dumpSuperKmers",
		"{p:jfPath} dump -c {i:jfDB} -o {o:kmerCount}")
	dumpSuperKmers.InParam("jfPath").FromStr(*jellyfishPath)
	dumpSuperKmers.SetOut("kmerCount", "{i:jfDB|%.jf}.txt")

	dumpSuperKmers.In("jfDB").From(countSuperKmers.Out("jfDB"))

	// construct the BFT
	constructBFT := wf.NewProc("bft",
		"{p:bftPath} build {p:k} kmers {i:colorFile} {o:bftOut}")
	constructBFT.InParam("bftPath").FromStr(*bftPath)
	constructBFT.InParam("k").FromInt(*kmerSize)
	constructBFT.SetOut("bftOut", filepath.Join(absGenomePath, "bft.{p:k}.out"))

	constructBFT.In("colorFile").From(colorFile.Out("colorFile"))

	// run kleuren
	runKleuren := wf.NewProc("kleuren",
		"{p:kleurenPath} -f {i:bftOut} -k {i:kmerCount} -b {o:bubblePath} -n {p:numMinColors} -d {p:maxDepth}")
	runKleuren.InParam("kleurenPath").FromStr(*kleurenPath)
	runKleuren.InParam("numMinColors").FromInt(*numMinColors)
	runKleuren.InParam("maxDepth").FromInt(*maxDepth)
	runKleuren.InParam("kmerSize").FromInt(*kmerSize)
	runKleuren.SetOut("bubblePath", filepath.Join(absGenomePath,
		"bubbles.kmer-{p:kmerSize}.depth-{p:maxDepth}."+time.Now().Format(time.Kitchen)+".out"))

	runKleuren.In("bftOut").From(constructBFT.Out("bftOut"))
	runKleuren.In("kmerCount").From(dumpSuperKmers.Out("kmerCount"))

	// Run workflow
	wf.Run()
}
