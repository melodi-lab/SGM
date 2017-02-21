#!/usr/bin/env python
#
# Copyright 2011 <fill in later>

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

"""Support for ionization, converting peptides to ions.

The second mass spectrometry step in MS/MS involves fragmenting a short
charged peptide, a precursor ion, into smaller fragments (product ions).
The most common product ions are b-ions and y-ions. We refer the reader
to [1] for a background on tandem mass spectrometry.

[1] Steen and Mann, The ABC's (and XYZ's) of Peptide Sequencing.
Available at https://j.ee.washington.edu/trac/msms/attachment/wiki/NewMembers/abc_of_peptide_seq.pdf

"""
from protein.peptide import Peptide, mass_table
from util.file import CommentedFile
from which.which import which

from subprocess import Popen, PIPE
import csv

class IonPredict(object):
    """Predict b/y ions using Crux command line tools.

    WARNING: This class uses the output of the program
    crux-predict-peptide-ions, which actually computes the ions given
    a peptide and precursor ion charge. If the output format changes,
    this class will no longer work. Tested against Crux v1.33.

    TODO(ajit): The tbl flag is currently ignored, and monoisotopic
    masses are used. Use the --fragment-mass flag on the command line
    to allow the user to choose which mass table gets used.

    Attributes:
        peptide: The peptide being fragmented.
        charge: The charge of the peptide as a precursor ion. If the charge
           is <= 0, then no ions will be returned.
        tbl: The ma

    Raises:
        WhichError: If crux-predict-peptide-ions is not in your $PATH.

    """
    def __init__(self, peptide, charge, tbl = 'mono'):
        """Predict peptide ions using crux-predict-peptide-ions.

        Arguments:
           peptide: Peptide sequence or instance of Peptide.
           charge: Charge of the peptide.
           tbl: Mass table to use, in (mono|average).

        """
        self.peptide = Peptide(str(peptide))
        self.charge = charge
        self.tbl = tbl

        cmd = [ 'crux-predict-peptide-ions' ]
        which(cmd[0]) # check that binary exists.
        args = [ str(self.peptide), '%d' % self.charge]
        self.reader = self._call(cmd + args)

    def _call(self, cmd):
        """Private function to run crux-predict-peptide-ions.

        Returns:
            A csv.DictReader which parses the non-comment lines output
            by crux-predict-peptide-ions.

        """
        f = Popen(cmd, stdout=PIPE, stderr=PIPE)
        cmts = [ '#', 'INFO:' ]
        self.ion_type = [ 'a', 'b', 'c', 'x', 'y', 'z' ]
        reader = csv.DictReader(CommentedFile(f.stdout, cmts), delimiter='\t')
        return reader

    @classmethod
    def _ion_name(cls, seq, charge, nh3, h2o):
        """Generate string name of an ion.

        Arguments:
           cls: Class
           seq: Peptide sequence [string]
           charge: Ion charge [int]
           nh3: Number of immonium losses [int]
           h2o: Number of water losses [int]

        """
        name = r'%s%s' % (seq, '+'*charge)
        if nh3:
            name = name + '(%sNH3)' % ('-' * nh3)
        if h2o:
            name = name + '(%sH2O)' % ('-' * h2o)
        return name

    def next(self):
        """Iterate over ions, returning a (name, mz) pair.

        This function is sensitive to changes in the format of
        crux-predict-peptide-ions output. This iterator requires
        that the file conforms to the following grammar

        file := comment | header(content)+ | info
        comment := '#'.*
        info := 'INFO:'.*

        header is a tab-separated list of fields: m/z, mass, charge,
        ion-series, peptide-bond-index, NH3, H2O, ISOTOPE, FLANK.
        content is a tab-separated list of values: float, float, int,
        int, int, int, int, int, int.

        comment and info lines are ignored.

        Returns:
           A triplet (m/z, name, kind) where m/z is the position of the
           ion in Thompsons [float], where name is a string description
           of the ion, and where kind is the ion type.

           kind: one of 'a', 'b', 'c', 'x', 'y', 'z'
           name: e.g., EAK+, EAK++, EAK+(--NH3), EAK++(--NH3)(-H2O)

        """
        for ion in self.reader:
            mz = float(ion['m/z'])
            charge = int(ion['charge'])

            seq = str(self.peptide)
            series = int(ion['ion-series'])
            index = int(ion['peptide-bond-index'])
            if series == 1:
                ion_seq = seq[0:index]
            elif series == 4:
                seq = seq[::-1]
                ion_seq = seq[0:index]
                ion_seq = ion_seq[::-1]
            else:
                raise ValueError('Bad call')

            name = IonPredict._ion_name(ion_seq, charge,
                              int(ion['NH3']), int(ion['H2O']))
            kind = self.ion_type[int(ion['ion-series'])]
            return (mz, name, kind)

        raise StopIteration

    def __iter__(self):
        return self
