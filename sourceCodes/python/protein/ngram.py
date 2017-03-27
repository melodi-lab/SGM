#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

"""
Computes the empirical transition probabilities in an n-gram model.  
Then outputs these learnt transitions to a dense-cpt file for use in
gmtk.  dense-cpt ordering that is defined in 
"""

from __future__ import with_statement

import os
from protein.peptide import Peptide, amino_acids_to_indices, amino_acids

__authors__ = [ 'Adam Gustafson <amg81@ee.washington.edu>' ]
    

def load_peptides(filename):
    """
    Load a file of peptide sequences into a list of peptides.

    Arguments:
        filename:   A file with peptide sequences, one per line.

    Returns:
        A list of Peptide sequences, one for each line in filename.
    """
    pep_list = []
    with open(filename) as f:
        for line in f:
            try:
                p = Peptide(line.strip())
                pep_list.append(p.seq)
            except ValueError:  # handles unknown "U" amino acids
                pass
    return pep_list


def prior(pep_list, alpha = 0):
    """
    Computes the prior probabilities of the amino acids in the provided list
    of peptides.
    
    Arguments:
        pep_list:   List of peptide sequences.
        alpha:      Dirichlet parameter used (default 0: no smoothing)
                    (See e.g. "Additive Smoothing" on Wikipedia)
    Returns:
        List of prior probabilities for the amino acids in 
        protein.peptides.amino_acids(), in the respective order. 
    """
    
    assert alpha >= 0, "Must have non-negative Dirichlet constant"    
    num_AA = len(amino_acids())
    indexer = amino_acids_to_indices()
    
    count = [0 for row in range(num_AA)]
    for p in pep_list:
        for aa in p:
            count[indexer[aa]] += 1
    N = sum(count)
    d = num_AA
    for i, val in enumerate(count):
        count[i] = (float(val) + alpha)/(N + alpha*d)

    return count


def transition_matrix(pep_list, n = 2, alpha = 0):
    """
    Computes the Markov transition probability matrix for the provided list
    of peptides.
    
    Arguments:
        pep_list:   List of peptide sequences.
        n:          Markov Order.  For example, a bigram corresponds to n = 2,  
                    a trigram corresponds to n = 3, etc.  Default is n = 2.
        alpha:      Dirichlet parameter used (default 0: no smoothing)
                    (See e.g. "Additive Smoothing" on Wikipedia)
    Returns:
        An array of size num_AA x (n x num_AA), where num_AA is the number of
        amino acids specified by protein.peptide.amino_acids() (currently 20).
        The (i,j)th entry of the matrix is equivalent to:
        P(AA_{t} = i | AA_{t-1} = i_1, AA_{t-2} = i_2, ...., AA_{t-m} = i_m),
        where j = 20^(m-1) i_1 + 20^(m-2) i_2 + .... + 20 i_{m-1} + i_m,
        where m = n - 1 is the length of the amino acids sequence in the 
        state-space.
    """
    # count:    will be transition matrix (initialize to zeros)
    #           (A list of n lists, each of length num_AA)
    #           The jth list is the jth column of the transition matrix
    
    assert alpha >= 0, "Must have non-negative Dirichlet constant"
    if n > 1:
        num_AA = len(amino_acids())
        count = [[0 for col in range(pow(num_AA,n-1))] for row in range(num_AA)]
    else:
        raise(ValueError, 'must have n > 1 (see ngram.prior for n = 1)')

    # iterator to get subsequences of length n in peptide sequence
    def subseq_iterator(sequence, n): 
        for i in range(0, len(sequence) - n + 1):
            yield list(sequence[i:(i + n)])
    
    # method to find the correct indices in the trans. matrix given the state
    def find_indices(obs, state, n):
        indexer = amino_acids_to_indices()
        m = n - 1
        i = indexer[obs]
        j = 0
        for k, ch in enumerate(state):
            j += pow(num_AA, (m - 1) - k) * indexer[ch]
        return i, j

    # count the occurrences of each state, and store in "count"
    for p in pep_list:
        for subseq in subseq_iterator(p, n):
            # Defining time to increase as move to the right... reverse subseq.
            # (interested in P(AA_t | AA_{t-1}, AA_{t-2}, ...) )
            subseq.reverse()            
            obs = subseq[0]
            state = subseq[1:]
            # print obs, state # for debugging
            i, j = find_indices(obs, state, n)
            count[i][j] = count[i][j] + 1

    # convert to transition probability matrix by dividing by row sums
    d = len(count[0])
    for i in range(num_AA):
        N = sum(count[i]) 
        if alpha > 0: # smoothing
            count[i] = [(v + alpha)/(N + alpha*d) for v in count[i]]
        elif N > 0: # no smoothing
            count[i] = [float(v)/N for v in count[i]]
        else: # no transitions to this amino acid -- will be row of zeros
            count[i] = [float(v) for v in count[i]]  
        
    return count    
    
    
def save_cpt(file_obj, cpt):
    """
    Writes a CPT to file object
    """
    
    def vec2sci_str(vec):
        sci_str = ''
        for val in vec:
            sci_str += '%e' % val + ' '
        sci_str += '\n'
        return sci_str
    
    # write matrix to file one line at a time 
    if isinstance(cpt[0], float): # have a vector, write the line
        file_obj.write(vec2sci_str(cpt))
    else:
        for row in cpt:
            save_cpt(file_obj, row)

            
def save_dense_cpt(filename, cpt_list, cpt_name_list, alpha_list):
    """
    Writes the Markov transition matrix to a CPT file for use in gmtk.

    Arguments:
        filename:           Where to write the cpt (e.g., 'dense-cpt.out')
        cpt_list:           A list of cpts calculated from ngram.prior
                            or ngram.transition_matrix
        cpt_name_list:      Name of the file in which to store cpt
        alpha_list:         List of dirichlet constants used 
    """
    if len(cpt_list) != len(cpt_name_list) or \
        len(cpt_list) != len(alpha_list) or \
        len(cpt_name_list) != len(alpha_list):
        raise(ValueError, 'All lists input should be of the same length.')
    
    n = len(cpt_list)
    num_AA = len(amino_acids())
    with open(filename, 'w') as f:
        f.write('% Dense CPTs \n')
        f.write(str(n) + ' % number of CPTs \n')
        card = ''
        for i in range(0, n):
            f.write(str(i) + ' % CPT #' + str(i) + ' of ' + str(n-1) + '\n')
            f.write(cpt_name_list[i] + '\n')
            f.write(str(i) + ' % number of parents \n')
            card += str(num_AA) + ' '
            f.write(card + ' % cardinalities \n')
            f.write('DirichletConst ' + str(alpha_list[i]) + '\n')
            save_cpt(f, cpt_list[i])

            
if __name__ == '__main__':
    dirpath = os.path.dirname(__file__)
    infile = dirpath + '/testdata/test.pep'
    outfile = dirpath + '/testdata/test-dense-cpt.out'
    pep_list = load_peptides(infile)
    n = 3
    cpt_list = []
    cpt_name_list = ['prologue_amino_acid_prior', 'padframe_amino_acid_prior',\
                       'amino_acid_prior']
    alpha_list = [0, 0, 0] # gmtk has [0.05, 0.0025, 0.000125]
    cpt_list.append(prior(pep_list, alpha_list[0]))
    for i in range(2, n+1): # i = 2 (bigram), i = 3 (trigram)
        cpt_list.append(transition_matrix(pep_list, i, alpha_list[i-1]))
    save_dense_cpt(outfile, cpt_list, cpt_name_list, alpha_list)