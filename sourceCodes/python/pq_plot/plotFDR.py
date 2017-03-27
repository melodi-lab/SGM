# From /homes/aklammer/riptide-dist/python/plotFDR.py

from pylab import *
import random

MAX_FDR = 0.10
random.seed(0)

#-------------------------------------------------------------------------------
def calcQ(positives, negatives, label="Notitle"):
  """Calculates P vs q xy points from a list of positives and negative scores"""
  all = []
  for positive in positives: all.append( (positive, 1) ) 
  for negative in negatives: all.append( (negative, 0) ) 

  random.shuffle(all)

  #--- sort descending
  all.sort( lambda x,y: cmp(y[0], x[0]) )
  ps = []
  fdrs = []
  posTot = 0.0
  fpTot = 0.0
  fdr = 0.0

  #--- iterates through scores
  for item in all:
    if item[1] == 1: posTot += 1.0
    else: fpTot += 1.0
    
    #--- check for zero positives
    if posTot == 0.0: fdr = 100.0
    else: fdr = fpTot / posTot

    #--- note the q
    fdrs.append(fdr)
    ps.append(posTot)

  qs = []
  lastQ = 100.0
  for idx in range(len(fdrs)-1, -1, -1):
    
    q = 0.0
    #--- q can never go up. 
    if lastQ < fdrs[idx]:
      q = lastQ
    else:
      q = fdrs[idx]
    lastQ = q
    qs.append(q)
  
  qs.reverse()
  # pos = range(len(qs))
  return qs, ps


#-------------------------------------------------------------------------------
def plotQfromFiles(positiveFile, negativeFile, color="g", style="-", \
    label="Notitle"):
  """ Plots P vs q from a two files of scores. To save a figure, use 
  e.g. pylab.savefig("Figure.eps") after calling this routine."""
  
  positives = map(float, open(positiveFile).readlines())
  negatives = map(float, open(negativeFile).readlines())

  plotQ(positives, negatives, color=color, style=style, label=label)

#-------------------------------------------------------------------------------
def plotQ(positives, negatives, color="g", style="-", label="Notitle"):
  """ Plots P vs q from a list of scores. To display a figure, use
  pylab.show(). To save a figure, use e.g. pylab.savefig("Figure.eps") after
  calling this routine."""

  qs, ps = calcQ(positives, negatives)

  newQs = []
  newPs = []
  for qIdx in range(len(qs)):
    if qs[qIdx] > MAX_FDR: 
      newQs.append(qs[qIdx])
      newPs.append(ps[qIdx])
      break
    newQs.append(qs[qIdx])
    newPs.append(ps[qIdx])

  for idx in range(len(newQs)):
    if newQs[idx] > 0.05:
      print >>sys.stderr, label + "\t" + str(newPs[idx])
      break
  plot(newQs,newPs, linestyle=style, color=color, label=label, linewidth=4)
  xlabel("q value", fontsize=24)
  ylabel("Positives", fontsize=24)
  xlim(0.0, MAX_FDR)

#-------------------------------------------------------------------------------
def calcBenjamini(scores):
  """Calculates the Benjamini-Hochberg FDR values and number of positives at
  that FDR from a list of p-values"""
  scores.sort()

  fdrs = [ i /1000.0 for i in range(1,100,1) ]
  positives = []
  for fdr in fdrs:
    # print >>sys.stderr, "At %.3f" % fdr
    jalphas = [ fdr * i / len(scores)  for i in range(1,len(scores)+1,1) ]
    for idx in range(len(scores)):
      if scores[idx] > jalphas[idx]:
        positives.append(idx-1)
        break
  return fdrs, positives

#-------------------------------------------------------------------------------
def plotBenjamini(positives, color="g", style="-", label="Notitle"):
  """Plots the Benjamini-Hochberg FDR values and number of positives for a
  list of p-values"""

  qs, ps = calcBenjamini(positives)


  for idx in range(len(qs)):
    if qs[idx] > 0.05:
      print >>sys.stderr, label + "\t" + str(ps[idx])
      break

  newQs = []
  newPs = []
  foundFDR = False
  for qIdx in range(len(qs)):
    if qs[qIdx] > MAX_FDR: 
      newQs.append(qs[qIdx])
      newPs.append(ps[qIdx])
      break
    if qs[qIdx] > 0.01 and foundFDR == False:
      positive = ps[qIdx] 
      foundFDR = True
    newQs.append(qs[qIdx])
    newPs.append(ps[qIdx])

  plot(newQs,newPs, linestyle=style, color=color, label=label)
  xlim(0.0, MAX_FDR)

#-------------------------------------------------------------------------------
def plotFDR(positives, negatives, label="Notitle", correctionFactor=1):
  all = []
  for positive in positives: all.append( (positive, 1) ) 
  for negative in negatives: all.append( (negative, 0) ) 

  all.sort( lambda x,y: cmp(y[0], x[0]) )
  fdrs = []
  ps = []
  posTot = 0
  fpTot = 0
  for item in all:
    if item[1] == 1: posTot += 1
    else: fpTot += (1 * correctionFactor)
    p = float(posTot)
    fdr = fpTot / p
    # give a little slack
    if fdr > MAX_FDR * 1.1: break
    fdrs.append(fdr)
    ps.append(p)
  plot(fdrs,ps,label=label)
  xlim(0, MAX_FDR)


if __name__ == '__main__':
  if len(sys.argv) < 4:
    raise SystemExit, "Usage: %s <pos-1> <neg-1> <pos-2> <neg-2>" % sys.argv[0]

  file1 = open(sys.argv[1])
  scores1 = []
  for line in file1:
    score = float(line.rstrip('\n'))
    scores1.append(score)
  file1.close()

  file2 = open(sys.argv[2])
  scores2 = []
  for line in file2:
    score = float(line.rstrip('\n'))
    scores2.append(score)
  file2.close()

  file3 = open(sys.argv[3])
  scores3 = []
  for line in file3:
    score = float(line.rstrip('\n'))
    scores3.append(score)
  file3.close()

  file4 = open(sys.argv[4])
  scores4 = []
  for line in file4:
    score = float(line.rstrip('\n'))
    scores4.append(score)
  file4.close()

  plotFDR(scores1, scores2, label="Sequest")
  plotFDR(scores3, scores4, label="Shark")
  legend()
  show()
