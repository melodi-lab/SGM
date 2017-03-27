#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

"""Peptide-spectrum matches, i.e., an assignment of a peptide to a spectrum.
"""

from __future__ import with_statement

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os

from protein.peptide import Peptide, amino_acids_to_indices

class PeptideSpectrumMatches(object):

    def __init__(self, psm_file = None, ignore_dup_spectra = False):
        """Initialize a peptide-spectrum match instance.

        Instances of this class do not explictly keep track of spectra and
        peptide objects. Rather they keep track of the pairing between spectrum
        ids, unique integer identifiers for spectra in a data set, and peptides.

        Args:
            psm_file: Name of file containing peptide <--> spectrum id matches.

        """
        self._spectrum_to_peptide = { }
        self._peptide_to_spectrum = { }
        if psm_file:
            self.parse_psm_file(psm_file, ignore_dup_spectra)

    def __len__(self):
        """Return the number of peptide-spectrum matches.

        """
        return len(self._spectrum_to_peptide)

    def __getitem__(self, key):
        """Retrieve a peptide-spectrum match, by peptide or spectrum-id.

        Args:
            key: Either a non-negative integer (spectrum id) or Peptide.

        Returns:
            If key is an integer, it returns the peptide match. If the key
            is a Peptide, it returns the spectrum match.

        Raises:
            TypeError: If the key is neither an 'int' nor a 'Peptide'.
            KeyError: If there is no match for the key.

        """
        typestring = type(key).__name__
        if typestring == 'int':
            return self._spectrum_to_peptide[key]
        elif typestring == 'Peptide':
            return self._peptide_to_spectrum[key]
        else:
            raise TypeError('Key has wrong type: %s' % typestring)

    def __setitem__(self, key, value):
        """Set a peptide-spectrum match, by peptide or spectrum.

        Treat peptide-spectrum matches as non-overwriteable. You can add
        a peptide-spectrum match, but you cannot overwrite or delete a
        peptide-spectrum match once it is added to the container. This is
        to prevent the user from shooting themselves in the foot, and ensures
        that there is a bijection between peptides and spectra.

        Args:
            key, value: One must be an integer, the other a Peptide.

        Raises:
            TypeError: If the key and value do not conform.
            KeyError: If you try to overwrite a peptide-spectrum match.

        """
        key_type = type(key).__name__
        value_type = type(value).__name__

        if key_type == 'int' and value_type == 'Peptide':
            if key in self._spectrum_to_peptide:
                raise KeyError('Spectrum %d already matched' % key)
            if value in self._peptide_to_spectrum:
                raise KeyError('Peptide %d already matches' % str(value))
            self._spectrum_to_peptide[key] = value
            self._peptide_to_spectrum[value] = key
        elif key_type == 'Peptide' and value_type == 'int':
            if key in self._peptide_to_spectrum:
                raise KeyError('Peptide %s already matched' % str(key))
            if value in self._spectrum_to_peptide:
                raise KeyError('Spectrum %d already matched' % value)
            self._peptide_to_spectrum[key] = value
            self._spectrum_to_peptide[value] = key
        else:
            raise TypeError('Invalid key and/or value type: %s %s' %
                            (key_type, value_type))

    def __iter__(self):
        return iter(self._spectrum_to_peptide.items())

    def peptides(self):
        """Return a list of peptides in the object."""
        return [ v for k, v in self ]

    def parse_psm_file(self, filename, ignore_dup_spectra):
        """Parse a psm-formatted file.

        PSM-formatted files store each peptide-spectrum match on a line.
        The first entry is an integer, the spectrum id. The second entry
        is the peptide match, written using IUPAC single letter codes.
        Entries are separated by whitespace characters.

        Args:
            filename: Name of the psm formatted file (e.g., foo.psm)

        Raises:
            RuntimeError: If any of the spectrum ids are negative, or
                if any of the lines contain more than two tokens.
            TypeError: If any of the spectrum ids are not integers.
            ValueError: If any of the peptides are not valid IUPAC coded
                peptide strings.

        """
        with open(filename) as psm_file:
            for lno, line in enumerate(psm_file):
                tokens = line.split()
                if len(tokens) != 2:
                    raise RuntimeError('Line %d too many entries: %s' %
                                       (lno, line))
                spectrum_id = int(tokens[0])
                peptide = Peptide(tokens[1])
                if spectrum_id < 0:
                    raise RuntimeError('Line %d bogus spectrum id: %s' %
                                       (lno, line))
                if ignore_dup_spectra:
                    if not spectrum_id in self._spectrum_to_peptide:
                        self[spectrum_id] = peptide
                else:
                    self[spectrum_id] = peptide

