/******************************************************************************
* test wrapper
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "wrap.h"

/******************************************************************************
* _hello()
* A test function for wrapping
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
*
******************************************************************************/
int wrhello( double **x, double *y, int nr, int nc)
{ 
    int j, k;
    double p;

    for (j = 0; j<nr; j++){
        for (k = 0; k< nc; k++){
            p = x[j][k] + y[k];
            printf("Val = %6.3f", p);
            printf(", Floor = %i", (int) floor(p));
            printf(", Ceil  = %i\n", (int) ceil(p));
        }
    }

    p = -123.93487;
    printf("p=%6.6f\n",p);
    k = (int) p;
    printf("p cast = %i\n",k);
    p = 123.93487;
    printf("p=%6.6f\n",p);
    k = (int) p;
    printf("p cast = %i\n",k);

    return (0);
}
/*****************************************************************************/

