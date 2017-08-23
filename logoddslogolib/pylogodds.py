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
#  Python version of log odds C code.
#
#  Authors: Stephen Altschul, Yi-Kuo Yu and Aleksandar Stojmirovic
#

import sys
from math import lgamma
import numpy as np
import logodds_data as data


num_dmixtures = len(data.dmixtures)


def _check_input_sequence(x, msg):

    try:
        y = np.array(x, dtype=np.float64)
    except ValueError:
        raise TypeError(msg)
    if len(y) != 4 and len(y) != 20:
        raise TypeError(msg)
    return y


def info1(counts, freqs):
    """
    Maximum-likelihood log-odds information content
    """


    msg_c = "counts must be a sequence of floats with size exactly 4 or 20"
    msg_p = "freqs must be a sequence of floats with size exactly 4 or 20"

    # alpha                     Size of alphabet - 1
    # info		        Information in bits
    # count			Number of observations
    # icount                    Integral portion of count
    # fraction		        Interpolation fraction
    # ptr			Pointer to data table
    # bound			Lower bound for asymptotics
    # A,B			Asymptotic constants

    c = _check_input_sequence(counts, msg_c)
    p = _check_input_sequence(freqs,  msg_p)

    if (len(c) != len(p)):
        raise TypeError("counts and freqs must have the same size")

    #  Set alternative data for proteins and DNA
    if (len(c) == 20):
        alpha = 19
        A = 27.14338
        B = 41.192
        bound = 185		#  Bound for fractional counts
        ptr = np.array(data.info1_proteindata, dtype=np.float64)
    else:
        alpha = 3
        A = 0.67425
        B = 3.052
        bound = 72		#  Bound for integral counts
        ptr = np.array(data.info1_dnadata, dtype=np.float64)

    #  Calculate total counts for column
    count = c.sum()

    #  Calculate maximum-likelihood per-letter log-odds score, in bits
    if count == 0:
        info = 0.0
    else:
        info = -count * np.log(count)
        ix = (c != 0.0)
        info += (c[ix] * np.log(c[ix]/p[ix])).sum()
        info /= count * np.log(2.0)

        #  Normalize log-odds score, using interpolation and data in table
        if count < bound:
            icount = int(np.floor(count))
            fraction = count - icount
            info -= (1.0-fraction) * ptr[icount-1] + fraction * ptr[icount]
        #  Normalize log-odds score, using asymptotic formula
        else:
            info -= (0.5*alpha*np.log(count)/log(2.0)-A+B/np.sqrt(count))/count

    # Return information content of column, in bits/letter
    return info * np.log(2.0)


def bildp(dmnumber, counts):
    """
    Protein-BILD information content
    """

    msg_c = "counts must be a sequence of floats with size exactly 20"
    c = _check_input_sequence(counts, msg_c)
    # check_input_sequence allows size 4 but bildp only covers proteins
    if len(c) != 20:
        raise TypeError(msg_c)

    #  weight -  Weight of Dirichlet component
    #  alpha  -  Dirichlet parameters
    weight, alpha = data.get_dirichlet_mixture_data(dmnumber)

    Alpha = alpha.sum(axis=1)                   #  Concentration parameters

    #  Base amino acid background frequencies on Dirichlet mixture
    aafq = (weight[:, np.newaxis] * alpha / Alpha[:, np.newaxis]).sum(axis=0)

    #  Calculate total counts for column
    count = c.sum()

    #  Calculate BILD score per letter
    if count == 0:
        info = 0.0
    else:
        for m in xrange(len(weight)):
            logprob = np.log(weight[m])+lgamma(Alpha[m])-lgamma(Alpha[m]+count)
            for i in xrange(20):
                logprob += lgamma(alpha[m][i]+c[i])-lgamma(alpha[m][i])
            if m == 0:
                info = logprob
            elif logprob > info:
                info = logprob + np.log(1.0 + np.exp(info-logprob))
            else:
                info += np.log(1.0 + np.exp(logprob-info))
        info -= (c * np.log(aafq)).sum()
        info /= count

    return info


def bildn(counts, freqs, alpha):
    """
    DNA-BILD information content
    """

    msg_c = "counts must be a sequence of floats with size exactly 4"
    msg_p = "freqs must be a sequence of floats with size exactly 4"

    c = _check_input_sequence(counts, msg_c)
    p = _check_input_sequence(freqs, msg_p)
    if (len(c) != len(p) or len(c) != 4):
        raise TypeError("counts and freqs must have size 4")

    #  Calculate total counts for column
    count = c.sum()

    # Calculate BILD score per letter
    if count == 0:
        info = 0.0
    else:
        info = lgamma(alpha) - lgamma(alpha + count)
        for i in xrange(4):
            info += lgamma(alpha*p[i]+c[i]) - lgamma(alpha*p[i])-c[i]*np.log(p[i])
        info /= count

    return info


def schneider(counts, freqs, flag):
    """
    Schneider's information content

    flag: no correct==0; correct==1
    """

    msg_c = "counts must be a sequence of floats with size exactly 4 or 20"
    msg_p = "freqs must be a sequence of floats with size exactly 4 or 20"
    c = _check_input_sequence(counts, msg_c)
    p = _check_input_sequence(freqs,  msg_p)
    if (len(c) != len(p)):
        raise TypeError("counts and freqs must have the same size")
    alpha = len(c)

    #  Calculate total counts for column
    count = c.sum()
    if count == 0:
        return 0.0

    #  Calculate entropy difference
    score = 0.0
    for i in xrange(alpha):
        if c[i] > 0.0:
            score += c[i] / count * np.log(c[i] / count)

    if not flag:
        score += np.log(alpha)
    else:

        #  Correct by interpolation between |_count_| and |_count_| + 1
        n = int(count)
        fraction = count - n

        #  Correct by mean entropy at floor
        sum = 0
        for i in xrange(alpha):
            diff = np.log(1/p[i] - 1)
            logprob = n * np.log(p[i])
            for j in xrange(n, 1, -1):
                sum += np.exp(logprob) * j * np.log(j)
                logprob += np.log(j / (n-j+1.0)) + diff

        score -= (1-fraction) * (sum/n - np.log(n))

        #  Correct by mean entropy at ceiling
        if fraction > 0.0:
            n += 1
            sum = 0
            for i in xrange(alpha):
                diff = np.log(1/p[i] -1)
                logprob = n * np.log(p[i])
                for j in xrange(n, 1, -1):
                    sum += np.exp(logprob) * j * np.log(j)
                    logprob += np.log(j / (n-j+1.0)) + diff

            score -= fraction * (sum/n - np.log(n))

    # Return column score, in nats
    return score


def _check_ords(aligns, msg):

    # aligns must be a list of lists
    # We check that they all have the same length, and that they are all
    # integers

    try:
        lengths = map(len, aligns)
        width = lengths[0]
        if any(map(lambda x: x != width, lengths)):
            raise TypeError(msg)
        if not all([all(map(lambda x: isinstance(x, int), row))
                    for row in aligns]):
            raise TypeError(msg)
        ords = np.array(aligns, dtype=np.int)
    except:
        raise TypeError(msg)
    return ords


def _effective(observed, freq):

    assert len(freq) == 20
    MAXIND = 400.0

    if observed == 20.0:
        return MAXIND

    high = 1.0
    while True:
        high *= 2
        distinct = 20.0 - np.exp(high * np.log(1 - freq)).sum()
        if distinct >= observed:
            break

    low = high / 2
    for iter in xrange(20):
        new = (low + high) / 2
        distinct = 20.0 - np.exp(new * np.log(1 - freq)).sum()
        if distinct < observed:
            low = new
        else:
            high = new

    if new > MAXIND:
        new=MAXIND
    return(new);


def weighcounts(aligns, freqs):
    """
    Weigh by effective number of observations
    """

    msg_ords = "ords must be a sequence of sequences of equal sizes"
    msg_p = "freqs must be a sequence of floats with size exactly 20"
    msg_error = "error! Numseen= %f, Colnm= %d\n"

    p = _check_input_sequence(freqs, msg_p)
    if (len(p) != 20):
        raise TypeError(msg_p)
    ords = _check_ords(aligns, msg_ords)
    numseq, width = ords.shape
    counts = np.zeros((width, 20), dtype=np.float64)

    seq_index = np.zeros(numseq, dtype=np.int)
    num = np.zeros(width, dtype=np.int)
    seen = np.zeros(20, dtype=np.int)

    #   Iterate on multiple alignment columns and amino acids
    for i in xrange(width):
        for j in xrange(20):
            #	Determine which and how many sequences contain amino acid number j.
            #	If <2 seqs. contain the a.a., then effect. no. of observs. is
            #   the count
            seqs = 0
            for k in xrange(numseq):
                if ords[k][i] ==j:
                    seq_index[seqs] = k
                    seqs += 1

            if seqs < 2:
                counts[i][j] = seqs
            else:
                #  Calculate number of distinct amino acids seen in filled
                #  columns other than the one under consideration
                num[:] = 0
                Numcol = 0
                for ii in xrange(width):
                    if ii == i:
                        continue

                    # Set flags for whether each amino acid has been seen to 0
                    numseen = 0
                    seen[:] = 0
                    for jj in xrange(seqs):
                        letter = ords[seq_index[jj]][ii]
                        #  Consider only columns without gaps
                        if letter >= 20:
                            break
                        #  Tally new letters for the column
                        elif seen[letter] == 0:
                            numseen += 1
                            seen[letter] = 1
                    else:
                        num[ii] = numseen
                        Numcol += 1

                # Use formula to estimate effective number of observations
                if Numcol == 0:
                    counts[i][j] = 1
                else:
                    Colnum = int((Numcol + 1.000001) / 2)
                    num *= -1
                    num.sort()
                    num *= -1
                    Numseen = np.float64(num[:Colnum].sum())

                    if Numseen > 0 and Colnum > 0:
                        counts[i][j] = _effective(Numseen/Colnum, p)
                        q = Numseen/Colnum
                    else:
                        raise RuntimeError(msg_error % (Numseen, Colnum))
                    if counts[i][j] > seqs:
                        counts[i][j] = seqs

    return counts
