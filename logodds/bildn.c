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

/*	DNA-BILD information content code				*/
/*	Version 1.0							*/
/*	May 17, 2013							*/
/*	Program by Stephen Altschul and Yi-Kuo Yu			*/

#include <math.h>
#include <stdio.h>
double bildn(double *c, double *p, double alpha);


double  bildn(double *c, double *p, double alpha)
       /*  c: Letter counts		*/
       /*  p: Background frequencies	*/
       /*    alpha Dirichlet parameter		*/
{
        double	count;			/*  Number of observations	*/
        double	info;			/*  Information in bits		*/
 	int	i;

/*  Calculate total counts for column					*/
 	for (count=i=0;i<4;++i) {
                count+=c[i];
                //fprintf(stderr,"i=%d, count=%lf, p=%lf\n",i,c[i],p[i]);
        }


/*  Calculate BILD score per letter					*/
        if (count == 0) info = 0.0;
        else {
                info=lgamma(alpha)-lgamma(alpha+count);
                //fprintf(stderr,"info = %lf count=%lf lgamm(2.0) =%lf lgamma(27.0)=%lf\n",info, count, lgamma(2.0), lgamma(27.0));
                for (i=0;i<4;++i)
                        info+=lgamma(alpha*p[i]+c[i])-lgamma(alpha*p[i])-c[i]*log(p[i]);

                info/= count*log(2.0);

        }
        return(info*log(2.0));
}
