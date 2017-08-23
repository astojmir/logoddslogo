#!/usr/bin/env python

# -------------------------------- LogOddsLogo --------------------------------
#
# Based on WebLogo source code that was modified and extended at the NCBI
#
#  All NCBI contributions are under Public Domain. See LICENSE.txt.
#

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  (1) Redistributions of source code must retain the above copyright notice,
#  this list of conditions and the following disclaimer.
#
#  (2) Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and or other materials provided with the distribution.
#
#  (3) Neither the name of the University of California, Lawrence Berkeley
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors
#  may be used to endorse or promote products derived from this software
#  without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.

# Replicates README.txt

"""
LogOddsLogo (http://www.ncbi.nlm.nih.gov/CBBresearch/Yu/logoddslogo/) is a tool
for creating sequence logos from biological sequence alignments. Building on
the WebLogo (http://weblogo.threeplusone.com) source code, LogOddsLogo uses
per-observation multiple-alignment log-odds scores as measures of information
content at each position of a sequence logo.

It can be run on the command line, as a standalone webserver, as a CGI webapp,
or as a python library.

Please consult the manual for installation instructions and more information:
htdocs/manual.html. The manual is also located at
http://www.ncbi.nlm.nih.gov/CBBresearch/Yu/logoddslogo/manual.html.

For help on the command line interface run
    ./logoddslogo --help

To build a simple logo run
    ./logoddslogo  < cap.fa > logo0.eps

To run as a standalone webserver at localhost:8080
    ./logoddslogo --serve

To create a logo in python code:
    >>> from logoddslogolib import *
    >>> fin = open('cap.fa')
    >>> seqs = read_seq_data(fin)
    >>> data = LogoData.from_seqs(seqs)
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options)
    >>> fout = open('cap.eps', 'w')
    >>> eps_formatter( data, format, fout)


-- Distribution and Modification --
This package is distributed under the new BSD Open Source License.
All NCBI contributions are under Public Domain.
Please see the LICENSE.txt file for details on copyright and licensing.
The LogOddsLogo source code can be downloaded from
http://code.google.com/p/weblogo/

LogOddsLogo requires Python 2.7, and the python array package 'numpy'
(http://www.scipy.org/scipylib/download.html)

Generating logos in PDF or bitmap graphics formats require that the ghostscript
program 'gs' be installed. Scalable Vector Graphics (SVG) format also requires
the program 'pdf2svg'.

"""

import sys
import math
from string import Template
import numpy as na

import corebio
from corebio.utils import ArgumentError
from corebio.utils import resource_string
from corebio.utils import isfloat
from corebio.seq import unambiguous_protein_alphabet

import weblogolib
from weblogolib import GhostscriptAPI
from weblogolib import classic
from weblogolib import std_color_schemes
from weblogolib import default_color_schemes
from weblogolib import std_units
from weblogolib import std_sizes
from weblogolib import std_alphabets
from weblogolib import std_percentCG
from weblogolib import pdf_formatter
from weblogolib import jpeg_formatter
from weblogolib import png_formatter
from weblogolib import png_print_formatter
from weblogolib import svg_formatter
from weblogolib import txt_formatter
from weblogolib import formatters
from weblogolib import base_distribution
from weblogolib import equiprobable_distribution
from weblogolib import read_seq_data

try:
    import _logodds as logodds
except ImportError:
    import pylogodds as logodds

DEFAULT_DMNUMBER_N = 1.0
DEFAULT_DMNUMBER_P = 2


# ------ META DATA ------

__all__ = [ 'LogoOptions',
            'description',
            '__version__',
            'LogoFormat',
            'LogoData',
            'GhostscriptAPI',
            'std_color_schemes',
            'default_color_schemes',
            'classic',
            'std_units',
            'std_sizes',
            'std_alphabets',
            'std_percentCG',
            'pdf_formatter',
            'jpeg_formatter',
            'png_formatter',
            'png_print_formatter',
            'txt_formatter',
            'eps_formatter',
            'formatters',
            'default_formatter',
            'base_distribution',
            'equiprobable_distribution',
            'read_seq_data',
            'color',
            'colorscheme',
            'logomath',
            ]

description  = "Create sequence logos from biological sequence alignments."

__version__ = '1.0.3'
release_description = "LogOddsLogo %s" % __version__



## Using the Robinson Robinson Frequencies
amino_acid_composition = dict(
A = 0.087135727479, R = 0.040894944605, N = 0.040418278548, D = 0.046870296325,
C = 0.033468612677, Q = 0.038269289735, E = 0.049525516559, G = 0.088606655336,
H = 0.033621241997, I = 0.036899088289, L = 0.085361634465, K = 0.080483022246,
M = 0.014743987313, F = 0.039767240243, P = 0.050677889818, S = 0.069597088795,
T = 0.058530491824, W = 0.010489421950, Y = 0.029922503029, V = 0.064717068767 )

aa_composition = [amino_acid_composition[_k]
                  for _k in unambiguous_protein_alphabet]
weblogolib.aa_composition = aa_composition

class LogoOptions(weblogolib.LogoOptions) :
    """ A container for all logo formatting options. Not all of these
    are directly accessible through the CLI or web interfaces.

    To display LogoOption defaults:
    >>> from logoddslib import *
    >>> LogoOptions()

    All physical lengths are measured in points. (72 points per inch, 28.3 points per cm)

    String attributes:
        o creator_text      -- Embedded as comment in figures.
        o logo_title
        o logo_label
        o unit_name         -- See std_units for options. (Default 'bits')
        o yaxis_label       -- Defaults to unit_name
        o xaxis_label
        o fineprint         -- Defaults to LogOddsLogo name and version

    Boolean attributes:
        o show_yaxis
        o show_xaxis
        o show_ends
        o show_fineprint
        o show_boxes        -- Draw boxes around stack characters (default: True)
        o debug             -- Draw extra graphics debugging information.
        o rotate_numbers    -- Draw xaxis numbers with vertical orientation?
        o scale_width       -- boolean, scale width of characters proportional to ungaps
        o pad_right         -- Make a single line logo the same width as multiline logos (default: False)

    Other attributes:
        o stacks_per_line

        o yaxis_tic_interval
        o yaxis_minor_tic_ratio
        o yaxis_scale
        o xaxis_tic_interval
        o number_interval

        o stack_neg  # a new variable to specify the span of negative scores
        o ovline_interval

        o shrink_fraction       -- Proportional shrinkage of characters if show_boxes is true.

        o errorbar_fraction
        o errorbar_width_fraction
        o errorbar_gray

        o resolution             -- Dots per inch (default: 96). Used for bitmapped output formats

        o default_color
        o color_scheme

        o stack_width           --
        o stack_aspect_ratio    -- Ratio of stack height to width (default: 5)

        o logo_margin           -- Default: 2 pts
        o stroke_width          -- Default: 0.5 pts
        o tic_length            -- Default: 5 pts
        o stack_margin          -- Default: 0.5 pts

        o small_fontsize        -- Small text font size in points
        o fontsize              -- Regular text font size in points
        o title_fontsize        -- Title text font size in points
        o number_fontsize       -- Font size for axis-numbers, in points.

        o text_font
        o logo_font
        o title_font

        o first_index
        o logo_start
        o logo_end

    """

    def __init__(self, **kwargs) :

        weblogolib.LogoOptions.__init__(self, **kwargs)

        # Changes from WebLogo
        self.score_method = "BILD"
        self.show_fineprint = False
        self.fineprint =  "LogOddsLogo" + __version__
        self.show_errorbars = False
        self.reverse_stacks = False
        self.ovline = False
        self.weighcounts = True

# End class LogoOptions


class LogoFormat(weblogolib.LogoFormat) :

    def __init__(self, data, options= None) :
        """ Create a new LogoFormat instance
        """

        # Making sure that changed, extended options are loaded
        assert options is not None
        assert isinstance(options, LogoOptions)

        weblogolib.LogoFormat.__init__(self, data, options)

        self.stack_height = self.stack_width * self.stack_aspect_ratio
        self.stack_neg = 0.0
        self.ovline_start, self.ovline_end = data.ovline_interval

        min_e = 10.0
        max_e = -10.0
        for i in range(data.length) :
            if data.entropy[i] < min_e :
                min_e = data.entropy[i]
            if data.entropy[i] > max_e :
                max_e = data.entropy[i]
        if min_e < -0.01 :
             self.stack_neg = math.ceil(self.stack_height *(-min_e)/max_e)

        conversion_factor = std_units[self.unit_name]
        if conversion_factor and max_e > math.log(len(self.alphabet)):
            self.yaxis_scale = max_e * conversion_factor

        min_stack_neg = (2.0 * self.stack_height *
                        self.yaxis_minor_tic_interval / self.yaxis_scale)
        if self.stack_neg < min_stack_neg:
            self.stack_neg = min_stack_neg

        self.line_height = (self.stack_height + self.stack_neg +
                            self.line_margin_top + self.line_margin_bottom )

        self.logo_height = int(2*self.logo_margin + self.title_height \
            + self.xaxis_label_height + self.line_height*self.lines_per_logo)

# End class LogoFormat



# ------ Logo Formaters ------

def eps_formatter(logodata, format, fout) :
    """ Generate a logo in Encapsulated Postscript (EPS)"""

    substitutions = {}
    from_format =[
        "creation_date",    "logo_width",           "logo_height",
        "lines_per_logo",   "line_width",           "line_height",
        "line_margin_right","line_margin_left",     "line_margin_bottom",
        "line_margin_top",  "title_height",         "xaxis_label_height",
        "creator_text",     "logo_title",           "logo_margin",
        "stroke_width",     "tic_length",
        "stacks_per_line",  "stack_margin",
        "yaxis_label",      "yaxis_tic_interval",   "yaxis_minor_tic_interval",
        "xaxis_label",      "xaxis_tic_interval",   "number_interval",
        "fineprint",        "shrink_fraction",      "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray",    "small_fontsize",       "fontsize",
        "title_fontsize",   "number_fontsize",      "text_font",
        "logo_font",        "title_font",
        "logo_label",       "yaxis_scale",          "end_type",
        "debug",            "show_title",           "show_xaxis",
        "show_xaxis_label", "show_yaxis",           "show_yaxis_label",
        "show_boxes",       "show_errorbars",       "show_fineprint",
        "rotate_numbers",   "show_ends",            "stack_height",
        "stack_width", "stack_neg", "ovline_start","ovline_end"
        ]

    for s in from_format :
        substitutions[s] = getattr(format,s)

    substitutions["shrink"] = str(format.show_boxes).lower()


    # --------- COLORS --------------
    def format_color(color):
        return  " ".join( ("[",str(color.red) , str(color.green),
            str(color.blue), "]"))

    substitutions["default_color"] = format_color(format.default_color)

    colors = []
    for group in format.color_scheme.groups :
        cf = format_color(group.color)
        for s in group.symbols :
            colors.append( "  ("+s+") " + cf )
    substitutions["color_dict"] = "\n".join(colors)


    ovline_start, ovline_end = format.ovline_start, format.ovline_end

    data = []

    # Unit conversion. 'None' for probability units
    conv_factor = std_units[format.unit_name]

    data.append("StartLine")

    seq_from = format.logo_start- format.first_index
    seq_to = format.logo_end - format.first_index +1

    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to) :
        logo_index = seq_index + format.first_index
        stack_index = seq_index - seq_from

        if stack_index!=0 and (stack_index % format.stacks_per_line) ==0 :
            data.append("")
            data.append("EndLine")
            data.append("StartLine")
            data.append("")

        data.append("(%s) StartStack" % format.annotate[seq_index] )

        if conv_factor:
            stack_height = logodata.entropy[seq_index] * std_units[format.unit_name]
        else :
            stack_height = 1.0 # Probability

        s = zip(logodata.counts[seq_index], logodata.alphabet)
        def mycmp( c1, c2 ) :
            # Sort by frequency. If equal frequency then reverse alphabetic
            if c1[0] == c2[0] : return cmp(c2[1], c1[1])
            return cmp(c1[0], c2[0])

        s.sort(mycmp)
        if format.reverse_stacks: s.reverse()

        C = float(sum(logodata.counts[seq_index]))
        if C > 0.0 :
            fraction_width = 1.0
            ## only if you scale the width, we don't do that
            ## if format.scale_width :
            ##     fraction_width = logodata.weight[seq_index]
            # print >>sys.stderr, fraction_width
            for c in s:
                data.append(" %f %f (%s) ShowSymbol" % (fraction_width, c[0]*stack_height/C, c[1]) )

        data.append("EndStack")
        if ovline_start <= seq_index <= ovline_end:
            data.append("DrawOverline")

        data.append("")

    data.append("EndLine")
    substitutions["logo_data"] = "\n".join(data)

    # Create and output logo
    template = resource_string( __name__, 'template.eps', __file__)
    logo = Template(template).substitute(substitutions)
    print >>fout, logo



# Monkey patch weblogolib to ensure that all other formatters call the
# correct eps_formatter
weblogolib.eps_formatter = eps_formatter

# Patch formatters - only eps_formatter is changed
formatters['eps'] = eps_formatter
default_formatter = eps_formatter

def _parse_explicit_composition(composition, alphabet, weight):

    explicit = composition[1: -1]
    explicit = explicit.replace(',',' ').replace("'", ' ').replace('"',' ')
    explicit = explicit.replace(':', ' ').split()

    if len(explicit) != len(alphabet)*2 :
        raise ValueError("Explicit prior does not match length of alphabet")
    prior = -na.ones(len(alphabet), na.float64)

    try :
        for r in range(len(explicit)/2) :
            letter = explicit[r*2]
            index = alphabet.ord(letter)
            value = float(explicit[r*2 +1])
            prior[index] = value
    except ValueError :
        raise ValueError("Cannot parse explicit composition")

    if any(prior==-1.) :
        raise ValueError("Explicit prior does not match alphabet")

    prior/= sum(prior)
    prior *= weight
    return prior


def parse_prior(composition, alphabet, weight=1.0) :
    """
    Parse a description of the expected monomer distribution of a nucleotide
    sequence. For protein sequences the prior is set implicitly for BILD
    scores and assumed either correspond to Robinson-Robinson frequencies (None
    - default) or to an explicit distribution (same as for nucleotides - see
    below).

    Valid compositions for nucleotides:

    - None  :                Use 'equiprobable'
    - 'equiprobable' :      All monomers have the same probability.
    - a percentage, e.g. '45%' or a fraction '0.45':
                            The fraction of CG bases for nucleotide alphabets
    - a species name, e.g. 'E. coli', 'H. sapiens' :
                            Use the average CG percentage for the specie's
                            genome.
    - An explicit distribution,  e.g. {'A':10, 'C':40, 'G':40, 'T':10}
    """

    if weight is None:
        weight = 1.0
    if weight < 0:
        raise ValueError("Weight cannot be negative")

    # Protein prior is Robinson-Robinson unless set explicitly
    if alphabet == unambiguous_protein_alphabet:
        if composition is None:
            prior = weight * na.array(aa_composition, na.float64)
        elif composition[0] == '{' and composition[-1] == '}':
            prior = _parse_explicit_composition(composition, alphabet, weight)
        else:
            prior = weight * na.array(aa_composition, na.float64)

    # Nucleotide prior depends on input
    else:
        if composition is None:
            composition = 'equiprobable'
        comp = composition.strip()

        if comp.lower() == 'equiprobable':
            prior = weight * equiprobable_distribution(len(alphabet))

        elif comp in std_percentCG :
            prior = weight * base_distribution(std_percentCG[comp])

        elif comp[-1] == '%' :
            prior = weight * base_distribution( float(comp[:-1]))

        elif isfloat(comp) :
            prior = weight * base_distribution( float(comp)*100. )

        elif composition[0] == '{' and composition[-1] == '}' :
            prior = _parse_explicit_composition(composition, alphabet, weight)

        else :
            raise ValueError("Unknown or malformed composition: %s"% composition)

    if len(prior) != len(alphabet) :
        raise ValueError(
            "The sequence alphabet and composition are incompatible.")

    return prior


class LogoData(weblogolib.LogoData) :
    """The data needed to generate a sequence logo.

    - alphabet
    - length
    - counts  -- An array of character counts
    - entropy -- The relative entropy of each column
    - ovline_interval -- overline start and ending point

     """

    def __init__(self, length=None, alphabet=None, counts=None, entropy=None,
                 weight=None, ovline_interval=None):
        """Creates a new LogoData object"""

        self.length = length
        self.alphabet = alphabet
        self.counts = counts
        self.entropy = entropy
        self.weight = weight
        self.ovline_interval = ovline_interval

    @classmethod
    def from_counts(cls, alphabet, counts, prior, score_method="BILD",
                    dmnumber=None, rawcounts=None, ovline=False):
        """Build a LogoData object from counts."""

        assert len(alphabet) == len(prior)
        assert len(alphabet) in (4, 20)

        seq_length, _ = counts.shape
        ent = na.zeros(seq_length, na.float64)

        if rawcounts is None:
            rawcounts = [0,]*seq_length
            for i in range (0, seq_length):
                rawcounts[i] = sum(na.array(counts[i], na.float64))

        if score_method == "NML":
            for i in range (0, seq_length) :
                alpha = na.array(counts[i], na.float64)
                ent[i] = logodds.info1(alpha, prior)

        elif score_method == "SU":
            for i in range (0, seq_length) :
                alpha = na.array(counts[i], na.float64)
                ent[i] = logodds.schneider(alpha, prior, 0)

        elif score_method == "SC":
            for i in range (0, seq_length) :
                alpha = na.array(counts[i], na.float64)
                ent[i] = logodds.schneider(alpha, prior, 1)

        elif score_method == "BILD":
            for i in range (0, seq_length) :
                alpha = na.array(counts[i], na.float64)
                if len(alphabet) == 4:
                    # dmnumber for nucleotides is a single floating point
                    # Dirichlet concentration parameter
                    if dmnumber is None:
                        dmnumber = DEFAULT_DMNUMBER_N
                    else:
                        dmnumber = float(dmnumber)
                        if dmnumber <= 0.0:
                            dmnumber = DEFAULT_DMNUMBER_N
                    ent[i] = logodds.bildn(alpha, prior, dmnumber)

                elif len(alphabet) == 20:
                    # dmnumber for proteins is an integer denoting a particular
                    # Dirichlet mixture
                    if dmnumber is None:
                        dmnumber = DEFAULT_DMNUMBER_P
                    else:
                        dmnumber = int(dmnumber)
                    ent[i] = logodds.bildp(dmnumber, alpha)

        else:
            raise ValueError("Unrecognised scoring method")


        weight = na.array( na.sum(counts,axis=1) , float)
        weight /= max(weight)

        # One-dimensional Smith-Waterman algorithm to determine the boundaries
        # for the overline
        if not ovline:
            ovline_interval = (-1, -1)
        else:
            S_loc = [max(0.0, rawcounts[0]*ent[0])]
            S_max = [0.0]
            i_start = [0]
            i_end = [0]

            for i in xrange(1, seq_length):
                cum_score = S_loc[-1] + rawcounts[i]*ent[i]
                if cum_score > S_max[-1]:
                    S_max[-1] = cum_score
                    i_end[-1] = i
                if cum_score > 0:
                    S_loc.append(cum_score)
                else:
                    S_loc.append(0.0)
                    S_max.append(0.0)
                    i_start.append(i+1)
                    i_end.append(0)

            # ovline_interval should be set properly within the loop
            #  This is just a sentinel
            ovline_interval = (None, None)
            F_max = -1
            for max_score, start_, end_ in zip(S_max, i_start, i_end):
                if max_score > F_max:
                    F_max = max_score
                    ovline_interval = (start_, end_)


        return cls(seq_length, alphabet, counts, ent, weight, ovline_interval)


    @classmethod
    def from_seqs(cls, seqs, prior, score_method="BILD", dmnumber=None,
                  weighcounts=False, ovline=False):
        """Build a LogoData object from a SeqList, a list of sequences."""

        # This should have been checked and reported in parse_prior(). Here we
        # just ensure that prior is well-formed
        assert len(seqs.alphabet) == len(prior)

        # check that at least one sequence of length at least 1 long
        if len(seqs)==0 or len(seqs[0]) ==0:
            raise ValueError("No sequence data found.")

        # Check sequence lengths
        err_fmt = "Sequence number %d differs in length from the previous sequences"
        seq_length = len(seqs[0])
        for i,s in enumerate(seqs) :
            if seq_length != len(s) :
                raise ArgumentError(err_fmt % (i+1), 'sequences')

        counts = seqs.profile()

        # Weight counts for proteins
        if weighcounts and len(seqs.alphabet) == 20:
            counts = na.array(logodds.weighcounts(seqs.ords(), prior))

        return cls.from_counts(seqs.alphabet, counts, prior, score_method,
                               dmnumber, None, ovline)
