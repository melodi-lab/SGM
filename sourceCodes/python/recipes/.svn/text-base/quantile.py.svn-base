#!/usr/bin/env bash
#
# http://adorio-research.org/download/Python/rquantile.py

from math import modf, floor

def quantile(x, q,  qtype = 7, issorted = False):
    """
    Args:
       x - input data
       q - quantile
       qtype - algorithm
       issorted- True if x already sorted.
       
    Compute quantiles from input array x given q.For median,
    specify q=0.5.
    
    References:
       http://reference.wolfram.com/mathematica/ref/Quantile.html
       http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile
    
    Author:
	Ernesto P.Adorio Ph.D.
	UP Extension Program in Pampanga, Clark Field.
    """
    if not issorted:
        y = sorted(x)
    else:
        y = x
    if not (1 <= qtype <= 9): 
       return None  # error!

    # Parameters for the Hyndman and Fan algorithm
    abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3

            (0,   0, 0, 1), # California linear interpolation, R type 4
            (0.5, 0, 0, 1), # hydrologists method, R type 5
            (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
            (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
            (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
            (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
           ]

    a, b, c, d = abcd[qtype-1]
    n = len(x)
    g, j = modf( a + (n+b) * q -1)
    if j < 0:
        return y[0]
    elif j > n:
        return y[n]

    j = int(floor(j))
    if g ==  0:
       return y[j]
    else:
       return y[j] + (y[j+1]- y[j])* (c + d * g)    

def Test():
    x = [11.4, 17.3, 21.3, 25.9, 40.1, 50.5, 60.0, 70.0, 75]

    for qtype in range(1,10):
        print qtype, quantile(x, 0.35, qtype)

if __name__ == "__main__":    
    Test()

"""
When the test code runs, it outputs

1 25.9
2 25.9
3 21.3
4 21.99
5 24.29
6 23.6
7 24.98
8 24.06
9 24.1175

This matches the output of the following R code:

> x <- c(11.4, 17.3, 21.3, 25.9, 40.1, 50.5, 60.0, 70.0, 75)
> for (i in seq(1,9)) { print(c( i, quantile(x, 0.35, type=i))) }
      35%
 1.0 25.9
      35%
 2.0 25.9
      35%
 3.0 21.3
        35%
 4.00 21.99
        35%
 5.00 24.29
      35%
 6.0 23.6
        35%
 7.00 24.98
        35%
 8.00 24.06
            35%
 9.0000 24.1175
>                                                
"""
