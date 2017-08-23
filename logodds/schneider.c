/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*/

/*	Schneider entropy-difference column-score code			*/
/*	Version 1.4							*/
/*	August 16, 2013							*/
/*	Program by Stephen Altschul					*/

#include <math.h>

double	schneider(int alpha, double *c, double *p, int flag)
		/* alpha:  Size of alphabet		*/
    		/*  *c:  Letter counts		*/
       		/*  *p: Background frequencies	*/
		/*  flag: no correct==0; correct==1	*/
{
	int	i,j,n;
	double	sum,diff,logprob;
	double	score;			/*  Entropy difference		*/
	double	count;			/*  Number of observations	*/
	double	fraction;		/*  Interpolation fraction	*/

/*  Calculate total counts for column					*/

	for (count=i=0;i<alpha;++i) count+=c[i];
	if (count==0) return(0.0);

/*  Calculate entropy difference					*/

	for (score=i=0;i<alpha;++i)
		if (c[i]) score+=c[i]/count*log(c[i]/count);
	if (!flag) score+=log((double) alpha);
	else {

/*  Correct by interpolation between |_count_| and |_count_| + 1	*/

		n=count;
		fraction=count-n;

/*  Correct by mean entropy at floor					*/

		for (sum=i=0;i<alpha;++i) {
			diff=log(1/p[i]-1);
			logprob=n*log(p[i]);
			for (j=n;j>1;--j) {
				sum+=exp(logprob)*j*log((double) j);
				logprob+=log(j/(n-j+1.0))+diff;
			}
		}
		score-= (1-fraction)*(sum/n-log((double) n));

/*  Correct by mean entropy at ceiling					*/

		if (fraction) {
			++n;
			for (sum=i=0;i<alpha;++i) {
				diff=log(1/p[i]-1);
				logprob=n*log(p[i]);
				for (j=n;j>1;--j) {
					sum+=exp(logprob)*j*log((double) j);
					logprob+=log(j/(n-j+1.0))+diff;
				}
			}
			score-= fraction*(sum/n-log((double) n));
		}
	}

/*  Return column score, in nats					*/

	return(score);
}
