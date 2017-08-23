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

/*	Amino acid effective count calculation code			*/
/*	Version 1.0							*/
/*	July 5, 2013							*/
/*	Program by Stephen Altschul and Yi-Kuo Yu    	                */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define	MAXIND	400.0			/* Max. no. of indep. observs.	*/

//void qsort(void *base, int nitems, int size, int (*compar)(const void *, const void*));

float	effective(double observed, double *freq)
{
	int	iter,j;
	double	distinct,low,high,new;

	if (observed==20.0) return(MAXIND);

	high=1.0;
	do {
		high*=2;
		distinct=20.0;
		for (j=0;j<20;++j) distinct-=exp(high*log(1-freq[j]));
	}
	while (distinct<observed);

	low=high/2;
	for (iter=0;iter<20;++iter) {
		new=(low+high)/2;
		distinct=20;
		for (j=0;j<20;++j) distinct-=exp(new*log(1-freq[j]));
		if (distinct<observed) low=new;
		else high=new;
	}
	if (new>MAXIND) new=MAXIND;
	return(new);
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)b - *(int*)a );
}

void  weighcounts(int **ords,double **counts, double *freq, int width, int numseq)
{
	int	i,j,k,ii,jj,I,m;
	int	seqs,letter,numseen,Numcol, Colnum;
	int	*seq_index, *num;
	int	seen[20];
	double	Numseen;

        seq_index = (int *)calloc(sizeof(int), numseq);
	num = (int *)calloc(sizeof(int), width);

/*	Iterate on multiple alignment columns and amino acids			*/

	for (i=0;i<width;++i) for (j=0;j<20;++j) {
                counts[i][j] = 0.0; /*initialize to zero count */
/*	Determine which and how many sequences contain amino acid number j.	*/
/*	If <2 seqs. contain the a.a., then effect. no. of observs. is the count	*/

		for (seqs=k=0;k<numseq;++k) if (ords[k][i]==j) seq_index[seqs++]=k;

		if (seqs<2) counts[i][j]=seqs;
		else {

/*	Calculate number of distinct amino acids seen in filled columns		*/
/*	other than the one under consideration					*/
                        for (ii=0;ii<width;ii++) num[ii]=0;
			for (Numseen=Numcol=ii=0;ii<width;++ii) if (ii!=i) {

                                        /*	Set flags for whether each amino acid has been seen to 0		*/

				for (numseen=I=0;I<20;++I) seen[I]=0;

				for (jj=0;jj<seqs;++jj) {
					letter = ords[seq_index[jj]][ii];

/*	Consider only columns without gaps					*/

					if (letter>=20) break;

/*	Tally new letters for the column					*/

					else if (seen[letter]==0) {
						++numseen;
						seen[letter]=1;
					}
				}

				if (jj==seqs) {
                                        num[ii]=numseen;
					//Numseen+=numseen;
					++Numcol;
				}
			}

/*	Use formula to estimate effective number of observations		*/

			if (Numcol==0) counts[i][j]=1;
			else {
                                Colnum = ((int) (Numcol+1.000001)/2);
                                qsort(num, width, sizeof(int), cmpfunc);
                                Numseen = 0.0;
                                for (m=0;m<Colnum;m++) Numseen += num[m];
                                if (Numseen > 0 && Colnum > 0)
                                        counts[i][j]=effective(Numseen/Colnum, freq);
                                else {
                                        printf("error! Numseen= %lf, Colnm= %d\n", Numseen, Colnum);
                                        exit(1);
                                }
                                if (counts[i][j]>seqs) counts[i][j]=seqs;
			}
		}
	}
        //free(seq_index);
        //free(num);
}

/*	Use binary search to find the effective number of indep. observations	*/
