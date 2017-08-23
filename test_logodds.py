#!/usr/bin/env python

#
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Tests for logodds library written at the NCBI
#

import unittest
from pkg_resources import resource_stream
import numpy as np

from corebio import seq_io
from corebio.seq import unambiguous_dna_alphabet
from corebio.seq import unambiguous_protein_alphabet

# logoddslogolib/_logodds.so must be compiled. Run Makefile in logodds/
from logoddslogolib import _logodds
from logoddslogolib import pylogodds
from logoddslogolib import std_alphabets
from logoddslogolib import aa_composition
from logoddslogolib import read_seq_data
from logoddslogolib import parse_prior
from logoddslogolib import equiprobable_distribution

PROTEIN_DATA = ['Prot.fa']
DNA_DATA = ['cap.fa']

PROTEIN_PRIOR = "{'A': 0.087135727479, 'C': 0.033468612677, 'E': 0.049525516559, 'D': 0.046870296325, 'G': 0.088606655336, 'F': 0.039767240243, 'I': 0.036899088289, 'H': 0.033621241997, 'K': 0.080483022246, 'M': 0.014743987313, 'L': 0.085361634465, 'N': 0.040418278548, 'Q': 0.038269289735, 'P': 0.050677889818, 'S': 0.069597088795, 'R': 0.040894944605, 'T': 0.058530491824, 'W': 0.01048942195, 'V': 0.064717068767, 'Y': 0.029922503029}"

def testdata_stream( name ):
    return resource_stream(__name__, 'test_logodds/data/'+name)


def load_seqs(fin, weight=None, alphabet='protein'):

    alphabet = std_alphabets[alphabet]
    seqs = read_seq_data(fin, seq_io.read, alphabet, ignore_lower_case=False)
    prior = np.array(parse_prior(None, seqs.alphabet, weight), np.float64)
    return seqs, prior


def get_counts(seqs, prior, logodds_mod=_logodds, weighcounts=False):

    counts = seqs.profile()
    seq_length, A = counts.shape
    ords = seqs.ords()
    if weighcounts and A == 20:
        counts = np.array(logodds_mod.weighcounts(ords, prior))
    return counts


def get_entropy(counts, prior, entropy_func):

    seq_length = counts.shape[0]
    ent = np.zeros(seq_length, np.float64)
    for i in range (0, seq_length) :
        alpha = np.array(counts[i], np.float64)
        ent[i] = entropy_func(alpha, prior)
    return ent


class test_weighcounts(unittest.TestCase) :

    def test_protein_counts(self):

        for test_file in PROTEIN_DATA:
            fin = testdata_stream(test_file)
            seqs, prior = load_seqs(fin, weight=None, alphabet='protein')
            counts1 = get_counts(seqs, prior, _logodds, True)
            counts2 = get_counts(seqs, prior, pylogodds, True)
            self.assertEqual(np.abs(counts1-counts2).sum(), 0.0)


def compare_all(data, efuncs):

    for seqs, prior in data:
        counts1 = get_counts(seqs, prior, _logodds, True)
        counts2 = get_counts(seqs, prior, pylogodds, True)
        for ef1, ef2 in efuncs:
            ent1 = get_entropy(counts1, prior, ef1)
            ent2 = get_entropy(counts2, prior, ef2)
            yield np.abs(ent1-ent2).sum()



class test_logodds(unittest.TestCase) :

    prot_data = [load_seqs(testdata_stream(test_file), None, 'protein')
                 for test_file in PROTEIN_DATA]
    dna_data = [load_seqs(testdata_stream(test_file), None, 'dna')
                for test_file in PROTEIN_DATA]

    def test_info1(self):

        data = self.prot_data + self.dna_data
        efuncs = [(_logodds.info1, pylogodds.info1)]
        for res in compare_all(data, efuncs):
            self.assertAlmostEqual(res, 0.0)

    def test_schneider0(self):

        ef1 = lambda alpha, prior: _logodds.schneider(alpha, prior, 0)
        ef2 = lambda alpha, prior: pylogodds.schneider(alpha, prior, 0)
        data = self.prot_data + self.dna_data
        efuncs = [(ef1, ef2)]
        for res in compare_all(data, efuncs):
            self.assertAlmostEqual(res, 0.0)

    def test_schneider1(self):

        ef1 = lambda alpha, prior: _logodds.schneider(alpha, prior, 1)
        ef2 = lambda alpha, prior: pylogodds.schneider(alpha, prior, 1)
        data = self.prot_data + self.dna_data
        efuncs = [(ef1, ef2)]
        for res in compare_all(data, efuncs):
            self.assertAlmostEqual(res, 0.0)

    def test_bildn(self):

        efs1 = [lambda alpha, prior, x=x: _logodds.bildn(alpha, prior, float(x))
                for x in range(1,6)]
        efs2 = [lambda alpha, prior, x=x: pylogodds.bildn(alpha, prior, float(x))
                for x in range(1,6)]
        for res in compare_all(self.dna_data, zip(efs1, efs2)):
            self.assertAlmostEqual(res, 0.0)

    def test_bildp(self):

        efs1 = [lambda alpha, prior, i=i: _logodds.bildp(i, alpha)
                for i in range(_logodds.num_dmixtures)]
        efs2 = [lambda alpha, prior, i=i: pylogodds.bildp(i, alpha)
                for i in range(pylogodds.num_dmixtures)]
        for res in compare_all(self.prot_data, zip(efs1, efs2)):
            self.assertAlmostEqual(res, 0.0)


# Original WebLogo test that fail in LogOddsLogo
class test_parse_prior(unittest.TestCase) :

    def test_explicit_protein(self) :
        self.assertTrue( all(1.0*np.array(aa_composition)  ==
            parse_prior(PROTEIN_PRIOR,  unambiguous_protein_alphabet ) ) )

    def test_weight(self) :
        self.assertTrue( all(1.0*equiprobable_distribution(4)  ==
            parse_prior( 'equiprobable',  unambiguous_dna_alphabet ) ) )
        self.assertTrue( all(123.123*equiprobable_distribution(4)  ==
            parse_prior( 'equiprobable',  unambiguous_dna_alphabet , 123.123) ) )

    def test_explicit(self) :
        s = "{'A':10, 'C':40, 'G':40, 'T':10}"
        p = np.array( (10, 40, 40,10), np.float64)*1.0 / 100.
        self.assertTrue( all(
            p == parse_prior( s,  unambiguous_dna_alphabet ) ) )



if __name__ == '__main__':
    unittest.main()
