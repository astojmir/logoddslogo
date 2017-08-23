#!/usr/bin/env python

# -------------------------------- LogOddsLogo --------------------------------
#
# Based on WebLogo source code that was modified and extended at the NCBI
#
#  All NCBI contributions are under Public Domain. See LICENSE.txt.
#

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006, The Regents of the University of California, through
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

import os
import sys
from string import Template
import cgi as cgilib
import cgitb; cgitb.enable()

from StringIO import StringIO
from weblogolib.color import *
from weblogolib.colorscheme import ColorScheme
from weblogolib.colorscheme import ColorGroup

import logoddslogolib
from corebio.utils import *

from weblogolib._cgi import resource_string
from weblogolib._cgi import mime_type
from weblogolib._cgi import extension
from weblogolib._cgi import alphabets
from weblogolib._cgi import color_schemes
from weblogolib._cgi import composition
from weblogolib._cgi import Field
from weblogolib._cgi import string_or_none
from weblogolib._cgi import truth
from weblogolib._cgi import int_or_none
from weblogolib._cgi import float_or_none


def _return_logo(logo, format, download):
    #
    #  RETURN LOGO OVER HTTP
    #
    #  Refactored as a separate function here for readability

    print "Content-Type:", mime_type[format]
    # Content-Disposition: inline       Open logo in browser window
    # Content-Disposition: attachment   Download logo
    if download:
        print 'Content-Disposition: attachment; ' \
            'filename="logo.%s"' % extension[format]
    else:
        print 'Content-Disposition: inline; ' \
            'filename="logo.%s"' % extension[format]

    # Separate header from data
    print

    # Finally, and at last, send the logo.
    print logo.getvalue()


score_methods_opts = ['BILD', 'NML', 'SU', 'SC']
score_methods = dict(zip(score_methods_opts, score_methods_opts))

def _convert_na_dmnumber(s):
    x = float(s)
    if not x > 0.0:
        raise ValueError()
    return x

def _convert_prot_dnmumber(s):
    x = int(s[-1])
    return x


def _get_controls(logooptions, proteins=True):

    # A list of form fields.
    # The default for checkbox values must be False (irrespective of
    # the default in logooptions) since a checked checkbox returns 'true'
    # but an unchecked checkbox returns nothing.

    if proteins:
        dm_field = Field('dmnumber',
                         'dm%d' % logoddslogolib.DEFAULT_DMNUMBER_P,
                         _convert_prot_dnmumber,
                         options=["dm%d" % i for i in range(logoddslogolib.logodds.num_dmixtures)],
                         errmsg="Invalid Dirichlet mixture index.")
        ab_field = Field('alphabet','alphabet_protein', alphabets.get,
                         options=['alphabet_protein'],
                         errmsg="Unknown sequence type.")
        default_color_scheme = 'color_hydrophobicity'

    else:
        dm_field = Field('dmnumber',
                         '%.2f' % logoddslogolib.DEFAULT_DMNUMBER_N,
                         _convert_na_dmnumber,
                         errmsg="Invalid Dirichlet mixture parameter.")
        ab_field = Field('alphabet','alphabet_dna', alphabets.get,
                         options=['alphabet_dna', 'alphabet_rna'],
                         errmsg="Unknown sequence type.")
        default_color_scheme = 'color_base_pairing'

    controls = [
        Field( 'score_method', 'BILD', score_methods.get,
               options=score_methods_opts,
               errmsg="Unknown scoring method."),
        Field( 'weighcounts', False, truth),
        Field( 'ovline', False, truth),
        dm_field,

        Field( 'sequences', ''),
        Field( 'format', 'png', logoddslogolib.formatters.get ,
            options=['png_print', 'png', 'jpeg', 'eps', 'pdf', 'svg', 'logodata'] ,
            errmsg="Unknown format option."),
        Field( 'stacks_per_line', logooptions.stacks_per_line , int,
            errmsg='Invalid number of stacks per line.'),
        Field( 'stack_width','medium', logoddslogolib.std_sizes.get,
            options=['small', 'medium', 'large'], errmsg='Invalid logo size.'),
        ab_field,

        Field( 'unit_name', 'bits',
            options=[ 'probability', 'bits', 'nats', 'kT', 'kJ/mol',
                        'kcal/mol']),
        Field( 'first_index', 1, int_or_none),
        Field( 'logo_start', '', int_or_none),
        Field( 'logo_end', '', int_or_none),
        Field( 'composition', 'comp_equiprobable', composition.get,
            options=['comp_equiprobable','comp_CG',
            'comp_Celegans','comp_Dmelanogaster','comp_Ecoli',
            'comp_Hsapiens','comp_Mmusculus','comp_Scerevisiae'],
            errmsg= "Illegal sequence composition."),
        Field( 'percentCG', '', float_or_none, errmsg="Invalid CG percentage."),
        Field( 'show_errorbars', False , truth),
        Field( 'logo_title', logooptions.logo_title ),
        Field( 'logo_label', logooptions.logo_label ),
        Field( 'show_xaxis', False, truth),
        Field( 'xaxis_label', logooptions.xaxis_label ),
        Field( 'show_yaxis', False, truth),
        Field( 'yaxis_label', logooptions.yaxis_label, string_or_none ),
        Field( 'yaxis_scale', logooptions.yaxis_scale , float_or_none,
            errmsg="The yaxis scale must be a positive number." ),
        Field( 'yaxis_tic_interval', logooptions.yaxis_tic_interval ,
                float_or_none),
        Field( 'show_ends', False, truth),
        Field( 'show_fineprint', False , truth),
        Field( 'color_scheme', default_color_scheme, color_schemes.get,
            options=color_schemes.keys() ,
            errmsg = 'Unknown color scheme'),
        Field( 'color0', ''),
        Field( 'symbols0', ''),
        Field( 'desc0', ''),
        Field( 'color1', ''),
        Field( 'symbols1', ''),
        Field( 'desc1', ''),
        Field( 'color2', ''),
        Field( 'symbols2', ''),
        Field( 'desc2', ''),
        Field( 'color3', ''),
        Field( 'symbols3', ''),
        Field( 'desc3', ''),
        Field( 'color4', ''),
        Field( 'symbols4', ''),
        Field( 'desc4', ''),
        Field( 'ignore_lower_case', False, truth),
        Field( 'scale_width', False, truth),
        ]

    return controls


def send_form(controls, template, errors=[]) :


    substitutions = {}
    substitutions["version"] = logoddslogolib.release_description
    # Bug fix. Not sure why this default substitution isn't added automatically like everything else
    substitutions['color_custom'] = ''
    for c in controls :
        if c.options :
            for opt in c.options :
                substitutions[opt.replace('/','_')] = ''
            substitutions[c.value.replace('/','_')] = 'selected'
        else :
            value = c.value
            if value == None : value = 'auto'
            if value=='true':
                substitutions[c.name] = 'checked'
            elif type(value)==bool :
                if value :
                    substitutions[c.name] = 'checked'
                else :
                    substitutions[c.name] = ''
            else :
                substitutions[c.name] = str(value)
        substitutions[c.name+'_err']  = ''
    substitutions['logo_range_err'] = ''

    # Disable graphics options if necessary auxiliary programs are not installed.
    try:
        command = find_command('gs')
    except EnvironmentError:
        try:
            command = find_command('gswin32c.exe')
        except EnvironmentError:
            substitutions['png_print'] = 'disabled="disabled"'
            substitutions['png'] = 'disabled="disabled"'
            substitutions['jpeg'] = 'disabled="disabled"'
            substitutions['pdf'] = 'disabled="disabled"'
            substitutions['svg'] = 'disabled="disabled"'
            substitutions['eps'] = 'selected="selected"'
    try:
        command = find_command('pdf2svg')
    except EnvironmentError:
        substitutions['svg'] = 'disabled="disabled"'


    if errors :
        print >>sys.stderr, errors
        error_message = []
        for e in errors :
            if type(e) is str :
                msg = e
            elif len(e)==2:
                substitutions[e[0]+"_err"] = "class='error'"
                msg = e[1]
            else :
                msg = e[0]


            error_message +=  "ERROR: "
            error_message +=  msg
            error_message += ' <br />'

        substitutions["error_message"] = ''.join(error_message)
    else :
        substitutions["error_message"] = ""

    html = Template(template).safe_substitute(substitutions)
    print "Content-Type: text/html\n\n"
    print html


def main(htdocs_directory, proteins) :

    logooptions = logoddslogolib.LogoOptions()
    controls = _get_controls(logooptions, proteins)

    if proteins:
        template = resource_string("proteins_template.html", htdocs_directory)
    else:
        template = resource_string("nucleicacid_template.html", htdocs_directory)

    form = {}
    for c in controls :
        form[c.name] = c


    form_values = cgilib.FieldStorage()

    # Send default form?
    if len(form_values) ==0 or form_values.has_key("cmd_reset"):
        # Load default truth values now.
        form['show_errorbars'].value = logooptions.show_errorbars
        form['show_xaxis'].value = logooptions.show_xaxis
        form['show_yaxis'].value = logooptions.show_yaxis
        form['show_ends'].value = logooptions.show_ends
        form['show_fineprint'].value = logooptions.show_fineprint
        form['scale_width'].value = logooptions.scale_width
        form['weighcounts'].value = logooptions.weighcounts

        send_form(controls, template)
        return

    # Get form content
    for c in controls :
        c.value = form_values.getfirst( c.name, c.default)

    options_from_form = ['format', 'stacks_per_line', 'stack_width',
        'alphabet', 'unit_name', 'first_index', 'logo_start','logo_end',
         'composition',
        'show_errorbars', 'logo_title', 'logo_label', 'show_xaxis',
        'xaxis_label',
        'show_yaxis', 'yaxis_label', 'yaxis_scale', 'yaxis_tic_interval',
        'show_ends', 'show_fineprint', 'scale_width',
        'weighcounts', 'ovline', 'score_method', 'dmnumber',
                         ]


    errors = []
    for optname in options_from_form :
        try:
            value =  form[optname].get_value()
            if value is not None:
                setattr(logooptions, optname, value)
        except ValueError, err :
            errors.append(err.args)


    # Construct custom color scheme
    custom = ColorScheme()
    for i in range(0,5) :
        color = form["color%d"%i].get_value()
        symbols = form["symbols%d"%i].get_value()
        desc = form["desc%d"%i].get_value()

        if color :
            try :
                custom.groups.append(ColorGroup(symbols, color, desc))
            except ValueError, e :
                errors.append( ('color%d'%i, "Invalid color: %s" % color) )

    if form["color_scheme"].value == 'color_custom' :
        logooptions.color_scheme =  custom
    else :
        try :
            logooptions.color_scheme = form["color_scheme"].get_value()
        except ValueError, err :
            errors.append(err.args)

    sequences = None

    # FIXME: Ugly fix: Must check that sequence_file key exists
    # FIXME: Sending malformed or missing form keys should not cause a crash
    # sequences_file = form["sequences_file"]
    if form_values.has_key("sequences_file") :
        sequences = form_values.getvalue("sequences_file")
        assert type(sequences) == str

    if not sequences or len(sequences)  ==0:
        sequences = form["sequences"].get_value()

    if not sequences or len(sequences)  ==0:
        errors.append( ("sequences", "Please enter a multiple-sequence alignment in the box above, or select a file to upload."))



    # If we have uncovered errors or we want the chance to edit the logo
    # ("cmd_edit" command from examples page) then we return the form now.
    # We do not proceed to the time consuming logo creation step unless
    # required by a 'create' or 'validate' command, and no errors have been
    # found yet.
    if form_values.has_key("cmd_edit") or errors :
        send_form(controls, template, errors)
        return

    # We write the logo into a local buffer so that we can catch and
    # handle any errors. Once the "Content-Type:" header has been sent
    # we can't send any useful feedback
    logo = StringIO()
    try:
        comp = form["composition"].get_value()
        percentCG = form["percentCG"].get_value()
        ignore_lower_case = form_values.has_key("ignore_lower_case")
        if comp=='percentCG': comp = str(percentCG/100)

        from corebio.matrix import Motif

        try:
            # Try reading data in transfac format first.
            motif = Motif.read_transfac(StringIO( sequences), alphabet=logooptions.alphabet)
            prior = logoddslogolib.parse_prior(comp, motif.alphabet)
            data = logoddslogolib.LogoData.from_counts(motif.alphabet,
                                                       motif,
                                                       prior,
                                                       logooptions.score_method,
                                                       logooptions.dmnumber,
                                                       None,
                                                       logooptions.ovline)

        except ValueError, motif_err :
            seqs = logoddslogolib.read_seq_data(StringIO( sequences),
                                        alphabet=logooptions.alphabet,
                                        ignore_lower_case=ignore_lower_case
                                        )
            prior = logoddslogolib.parse_prior(comp, seqs.alphabet)
            data = logoddslogolib.LogoData.from_seqs(seqs, prior,
                                                     logooptions.score_method,
                                                     logooptions.dmnumber,
                                                     logooptions.weighcounts,
                                                     logooptions.ovline)

        logoformat =  logoddslogolib.LogoFormat(data, logooptions)
        format = form["format"].value
        logoddslogolib.formatters[format](data, logoformat, logo)
    except ValueError, err :
        errors.append( err.args )
    except IOError, err :
        errors.append( err.args)
    except RuntimeError, err :
        errors.append( err.args )

    if form_values.has_key("cmd_validate") or errors :
        send_form(controls, template, errors)
        return

    _return_logo(logo, format, form_values.has_key("download"))
