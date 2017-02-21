#!/usr/bin/python
# AUTHOR: William Stafford Noble
# CREATE DATE: May 3, 2010
#
# ORIGIN: http://noble.gs.washington.edu/~wnoble/proj/pivdo/results/bill/2010-05-03/make-pq-plot.py

import sys
import os

USAGE = """USAGE: make-pq-plot.py [options] <targets> <decoys> <root>

  Uses gnuplot to produce a png plot of the number of targets as a
  function of q-value.  The q-value is the minimum FDR at which a
  given score is deemed significant.  FDR of threshold x is estimated
  as number of decoys with score > x divided by number of targets with
  score > x.  If the number of targets and decoys is not equal, a
  normalizing factor is included in the FDR estimate.

  The following files are produced:

    <root>.txt - Five columns (#targets, score, #decoys, FDR, q-value)
    <root>.gnuplot - Gnuplot file for producing plot
    <root>.png - The plot itself

  OPTIONS:
    --column <int>   Column in which scores reside (indexed from 0)
    --fdr            Include FDR as a second series
    --header <int>   Number of header lines to ignore in the data files.

  Assumes that the input files are tab-delimited.

  If one of the input filenames is specified as "-", then stdin is read.
"""

##############################################################################
# Read the values for a specified input column into an array.
def readScoreFile(scoreColumn, scoreFileName, nHeaderLines):
  if (scoreFileName == "-"):
    scoreFile = sys.stdin
  else:
    scoreFile = open(scoreFileName, "r")

  returnValue = []
  lineNum = 0
  for line in scoreFile:
    line = line.rstrip()
    words = line.split()
    if (len(words) < scoreColumn):
      sys.stderr.write("Error parsing line %d.\n%s\n" % (lineNum, line))
      sys.exit(1)
    if (lineNum >= nHeaderLines):
      score = float(words[scoreColumn])
      returnValue.append(score)
    lineNum += 1
  scoreFile.close()

  sys.stderr.write("Read %d scores from %s.\n" %
                   (len(returnValue), scoreFileName))
  return(returnValue)
    

##############################################################################
# MAIN
##############################################################################
# Set default options.
scoreColumn = 0
includeFDR = 0
nHeaderLines=0

# Parse the command line.
sys.argv = sys.argv[1:]
while (len(sys.argv) > 3):
  nextArg = sys.argv[0]
  sys.argv = sys.argv[1:]
  if (nextArg == "--column"):
    scoreColumn = int(sys.argv[0])
    sys.argv = sys.argv[1:]
  elif (nextArg == "--fdr"):
    includeFDR = 1
  elif (nextArg == "--header"):
    nHeaderLines = int(sys.argv[0])
    sys.argv = sys.argv[1:]
  else:
    sys.stderr.write("Invalid option (%s).\n" % nextArg)
    sys.exit(1)
if (len(sys.argv) != 3):
  sys.stderr.write(USAGE)
  sys.exit(1)
targetFileName = sys.argv[0]
decoyFileName = sys.argv[1]
outputRoot = sys.argv[2]

# Read all the scores from the files.
targets = readScoreFile(scoreColumn, targetFileName, nHeaderLines)
decoys = readScoreFile(scoreColumn, decoyFileName, nHeaderLines)

# Sort the scores.
targets.sort(reverse=True)
decoys.sort(reverse=True)

# Factor to account for relative size of targets and decoys.
ratio = float(len(targets)) / float(len(decoys))
sys.stderr.write("Ratio = %g\n" % ratio)

# Compute FDRs.
numDecoys = []
fdrs = []
decoyIndex = 0
decoy = decoys[decoyIndex]
for targetIndex in xrange(0, len(targets)):
  target = targets[targetIndex]

  # Find the first decoy < this target score.
  while (decoy >= target):
    decoyIndex = decoyIndex + 1
    if (decoyIndex >= len(decoys)):
      break
    decoy = decoys[decoyIndex]

  # Estimate the FDR, with a ceiling of 1.0.
  fdr = ratio * (float(decoyIndex) / float(targetIndex + 1))
  if (fdr > 1.0):
    fdr = 1.0

  numDecoys.append(decoyIndex)
  fdrs.append(fdr)

# Compute q-values.
qvalues = []
maxFDR = 0.0
for fdrIndex in xrange(0, len(fdrs)):
  if (fdrs[fdrIndex] > maxFDR):
    maxFDR = fdrs[fdrIndex]
  qvalues.append(maxFDR)

# Print the data file.
outFileName = "%s.txt" % outputRoot
outFile = open(outFileName, "w")
for index in xrange(0, len(qvalues)):
  outFile.write("%d\t%g\t%d\t%g\t%g\n" 
                % (index, targets[index], numDecoys[index], 
                   fdrs[index], qvalues[index]))
outFile.close()

# Print the gnuplot file.
gnuplotFileName = "%s.gnuplot" % outputRoot
gnuplotFile = open(gnuplotFileName, "w")
gnuplotFile.write("set output '/dev/null'\n")
gnuplotFile.write("set xrange [0:0.1]\n")
gnuplotFile.write("set terminal png\n")
gnuplotFile.write("set xlabel 'q-value'\n")
gnuplotFile.write("set ylabel 'Number of decoys'\n")
if (includeFDR == 1):
  gnuplotFile.write("plot '%s.txt' using 4:1 notitle with lines\nre" % 
                    outputRoot)
gnuplotFile.write("plot '%s.txt' using 5:1 notitle with lines\n" %
                  outputRoot)
gnuplotFile.write("set output\n")
gnuplotFile.write("replot\n")
gnuplotFile.close()

# Create the plot.
command = "gnuplot %s.gnuplot > %s.png" % (outputRoot, outputRoot)
returnValue = os.system(command)
if (returnValue != 0):
  sys.stderr.write("Error running gnuplot.\n%s\n" % command)
  sys.exit(1)
