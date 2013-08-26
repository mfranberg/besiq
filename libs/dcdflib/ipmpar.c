#include <float.h>
#include <limits.h>

int ipmpar(int*);
/*
-----------------------------------------------------------------------
 
     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
 
  INTEGERS.
 
     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
 
               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
 
               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
 
     IPMPAR(1) = A, THE BASE.
 
     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
 
     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
 
  FLOATING-POINT NUMBERS.
 
     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
 
               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
 
               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
 
     IPMPAR(4) = B, THE BASE.
 
  SINGLE-PRECISION
 
     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
 
     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
 
     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
 
  DOUBLE-PRECISION
 
     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
 
     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
 
     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
 
-----------------------------------------------------------------------
 
     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED REMOVE
     THE COMMENT DELIMITORS FROM THE DEFINITIONS DIRECTLY BELOW THE NAME
     OF THE MACHINE
 
-----------------------------------------------------------------------
 
     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
 
-----------------------------------------------------------------------
     .. Scalar Arguments ..
*/
int ipmpar(int *i)
{
    static int imach[11] = {
        0, /* dummy zeroth index */
        /* integers */  
        2, /* base;  you are on the binary machine, aren't you? */
        sizeof(int)*(CHAR_BIT)-1,
        INT_MAX,
        /* all floats */
        FLT_RADIX,
        /* single precision floats */
        FLT_MANT_DIG,
        FLT_MIN_EXP,
        FLT_MAX_EXP,
        /* double precision floats */
        DBL_MANT_DIG,
        DBL_MIN_EXP,
        DBL_MAX_EXP,
    };
    return imach[*i];
}
