#include <stdio.h>
#include <math.h>
#include "gauss.h"

//Global matrix
float matrix[DEGREES][DEGREES];

/*<f>*/
/*@f*****************************************************************************
 *
 *  CONFIGURATION
 *      div_tab
 *
 *  DESCRIPTION
 *      Configuration table of divisor
 *
 *      DIV_TAB_STRUCT div_tab
 *
 *@f****************************************************************************/

DIV_TAB_STRUCT div_tab = {

/* threshold, default, warning channel, enable warning  */

   0.001,     0.001,     //W_FCS_DIV_16, &calc_config_data.warning_check_enable,   /* matrix[][] */

};

int main () {
	float x[7] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
	float y[7] = {1.0, 2.0, 0.0, 2.0, 3.0, 2.0, 8.0};
	float coef[DEGREES];
	int result = gauss(x, y, DEGREES, sizeof(x)/sizeof(x[0]), coef);
	for ( int i = 0; i < DEGREES; i++ )
	{
		printf("%lf\n", coef[i]);
	}

	return 0;
}

/*<f>*/
/*@f***************************************************************************
 *
 *  NAME
 *      gauss
 *
 *  SYNTAX
 *      int gauss (float *fx_ptr, float *fy_ptr, int no_st, float *coeff_ptr)
 *
 *      where:
 *          fx_ptr    : points to the x-positions of the valid stress_zones;
 *          fy_ptr    : points to the strip stress regarding to the valid
 *                      stress zones;
 *          no_st     : number of the active stress zones;
 *          coeff_ptr : points to the coefficients;
 *
 *  DESCRIPTION
 *          According to the x[] and y[] values the coefficients of this
 *          function are calculated by using the Gauss polynom approximation
 *
 *          INPUT
 *          
 *          fx[]              : x-positions of the valid stress zones between
 *                              [ -1  +1 ];
 *          fy[]              : strip stress regarding on the valid stress zones
 *                              in [kN/mm^2];
 *          no_st             : number of the active stress zones; 
 *
 *          OUTPUT
 *
 *          coeff.k.vector[0]  : coefficient a0;
 *          coeff.k.vector[1]  : coefficient a1;
 *          coeff.k.vector[2]  : coefficient a2;
 *          coeff.k.vector[3]  : coefficient a3;
 *          coeff.k.vector[4]  : coefficient a4;
 *
 *
 *  RETURNS
 *      this function returns the error code
 *
 *      error = 0 => poly_approx passed
 *              1 => number of active stress zones to small
 *              2 => no solution found
 *              3 => matrix is not symmetrical
 *              4 => matrix is singular
 *              5 => matrix is not defined positiv
 *              6 => mode is false
 *
 *@f**************************************************************************/

int gauss (float *fx_ptr, float *fy_ptr, int degrees, int no_st, float *coeff_ptr)
{

    float y[degrees];
    float x[degrees];
    int i, j, k;
    int ii,jj;
    float sh, sm;
    float s;
    float y_0;
    int error, ready;

    /* 
     * Initialization of x[]
     */
    for (i = 0;i < degrees;i++) 
            x[i] = 0.0;

    /* 
     * number of base points to small  
     */
     
    if (no_st < degrees)
        return(1);

    /* 
     * reset matrix and coefficient vector   
     */

    for (i = 0;i < degrees;i++) {
        for (j = 0;j < degrees;j++) {
            matrix[i][j] = 0.0;
        }
        coeff_ptr[i] = 0.0;
    }

    /* 
     * calculate upper triangle of the matrix   
     */

    for (k = 0;k < no_st;k++) {
        funk(fx_ptr[k],&y[0]);
        for (i = 0;i < degrees;i++) {

            for (j = 0;j < degrees;j++) {
                matrix[i][j] = matrix[i][j] + y[i] * y[j];
            }
            coeff_ptr[i] = coeff_ptr[i] + y[i] * fy_ptr[k];
        }
    }

    /* 
     * supplement the matrix
     */

    for (i = 1;i < (degrees + 1);i++) {
        ii = i - 1;
        for (j = 1;j < i;j++) {
            jj = j - 1;
            matrix[ii][jj] = matrix[jj][ii];
        }
    }

    /* 
     * solution of the equation
     */

    error = pos_def_matrix(coeff_ptr, degrees, (int)1);
    if (error != 0)
        return(error);
     

    /*
     * determination of the mean square factor before first iteration
     */

    sh = 0.0;
    for (i = 0;i < no_st;i++) {
        sh = sh + fy_ptr[i] * fy_ptr[i];
    }
    sm = sh;


    /* 
     * determination of a new mean square factor 
     */

    ready = FALSE;
    while(ready == FALSE) {
        s = 0.0;
        for (i = 0;i < degrees;i++) {
            x[i] = 0.0;
        }
        for (k = 0;k < no_st;k++) {
            funk(fx_ptr[k],&y[0]);
            y_0 = fy_ptr[k];
            for (i = 0;i < degrees;i++)
                y_0 = y_0 - coeff_ptr[i] * y[i];
            s = s + y_0 * y_0;
            for (i = 0;i < degrees;i++)
                x[i] = x[i] + y[i] * y_0;
        }
        if (s >= sm)
            ready = TRUE;
        if (ready == FALSE) {
            sm = s;
            error = pos_def_matrix(&x[0], degrees, (int)3);
            for (i = 0;i < degrees;i++)
                coeff_ptr[i] = coeff_ptr[i] + x[i];
        }
    }

    /*
     * if mean square error == sum of squares, then no solution found
     */

    if (fabs(sm - sh) < THRESHOLD)
        return(2);
    else
        return(0);
}



/*<f>*/
/*@f***************************************************************************
 *
 *  NAME
 *      funk
 *
 *  SYNTAX
 *      void funk (float x,float *y_ptr)
 *
 *      where:
 *          x     : the actual x position between [-1  +1];
 *          y_ptr : y-values;
 *
 *  DESCRIPTION
 *          Calculates the y-values depending on the actual x-position;
 *
 *          INPUT
 *
 *          x     : the actual x position between [-1  +1];
 *
 *          OUTPUT
 *
 *          y[]   : y - values depending on the actual x position;
 *
 *  RETURNS
 *      ----
 *
 *@f**************************************************************************/


void funk (float x,float *y_ptr)
{
    y_ptr[0] = 1;
    y_ptr[1] = x;
    y_ptr[2] = x * x;
    y_ptr[3] = x * x * x;
    y_ptr[4] = y_ptr[2] * y_ptr[2];

}


/*<f>*/
/*@f***************************************************************************
 *
 *  NAME
 *      pos_def_matrix
 *
 *  SYNTAX
 *      int pos_def_matrix (float *x_ptr,int n ,int mode)
 *
 *      where:
 *
 *  DESCRIPTION
 *
 *  RETURNS
 *      this function returns the error code
 *
 *      error = 0 => calculation passed
 *              3 => matrix is not symmetrical
 *              4 => matrix is singular
 *              5 => matrix is not defined positiv
 *              6 => mode is false
 *
 *@f**************************************************************************/


int pos_def_matrix ( float *x_ptr, int n, int mode)
{

    float s;
    int i,ii,j,k,kk;

    /*
     * only mode 1 u 3
     */

    if ((mode == 1) || (mode == 3)) {

        if (mode == 1) {

            for (i = 0;i < n;i++) {

                for (j = (i + 1);j < n;j++) {
                    if ((fabs(matrix[i][j]) - fabs(matrix[j][i])) > THRESHOLD)
                        return(5);
                }
            }
            for (i = 1;i < (n+1);i++) {
                ii = i - 1;
                s = matrix[ii][ii];
                for (k = 1;k < i;k++)
                    s = s - matrix[k-1][ii] * matrix[k -1][ii];
                if (fabs(s) < THRESHOLD)
                    return(3);
                if (s < 0.0)
                    return(4);
                matrix[ii][ii] = (float)(sqrt((double)s));
                for (j = i;j < n ;j++) {
                    s = matrix[ii][j];
                    for (k = 1;k < i;k++) {
                        kk = k - 1;
                        s = s - matrix[kk][ii] * matrix[kk][j];
                    }
                    matrix[ii][ii] = check_divisor(matrix[ii][ii],
                                     &div_tab.divisor[0]);
                    matrix[ii][j] = s / matrix[ii][ii];
                }
            }
        }
        for (i = 1;i < (n+1);i++) {
            ii = i - 1;
            s = x_ptr[ii];
            for (k = 1;k < i;k++) {
                kk = k - 1;
                s = s - matrix[kk][ii] * x_ptr[kk];
            }
            matrix[ii][ii] = check_divisor(matrix[ii][ii],
                                           &div_tab.divisor[0]);
            x_ptr[ii] = s / matrix[ii][ii];
        }
        for (i = 0;i < n;i++) {
            j = n - 1 - i;
            s = x_ptr[j];
            for (k = (j + 1);k < n ;k++)
                s = s - matrix[j][k] * x_ptr[k];
            matrix[j][j] = check_divisor(matrix[j][j],
                                         &div_tab.divisor[0]);
            x_ptr[j] = s / matrix[j][j];
        }
        return(0);
    }
    else
        return(6);
}

/*<f>*/
/*@f***************************************************************************
 *
 *  NAME
 *      check_divisor
 *
 *  SYNTAX
 *      float check_divisor(float divisor, DIV_STRUCT *div_ptr)
 *
 *      where:
 *          divisor - is the value which has to be checked  
 *          div_ptr - is a pointer to a DIV_STRUCT object  
 *
 *  DESCRIPTION
 *      If the absolute value of the divisor exceeds the threshold value of the
 *      div_tab (initialized in dbfields.c) this procedure returns the divisor.
 *      In the other case, the procedure returns the default value of the div_tab.
 *
 *  RETURNS
 *      divisor - if absolute value of divisor > threshold value of div_tab
 *      default - if absolute value of divisor < threshold value of div_tab
 *
 *@f**************************************************************************/


float check_divisor(float divisor, DIV_STRUCT *div_ptr)
{
/*    trace_msg((BYTE)0x5,(BYTE)5, "warning_check_enable: %d", *div_ptr->warning_check_enable);*/
    //if (*div_ptr->warning_check_enable) 
    //    alr_signal_alarm_channel(div_ptr->warning_channel, 
    //        (int)(fabs(divisor) < div_ptr->threshold));

    if (fabs(divisor) < div_ptr->threshold) 
        return(div_ptr->default_value);

    else
        return(divisor);
}
