/*@m****************************************************************************
 * ctrllib.h - define/typedefs and function prototypes for controller kernel
 * $Revision:   1.0  $
 * $Logfile:   C:/TSSC/MILLSTND/APPS/PROCON/FCS/VCS/CTRLLIB.H_V  $
 *@m****************************************************************************/


#ifndef _ctrllibh
#define _ctrllibh


/*
 *  prototypes
 */

#define THRESHOLD (float) 0.001
#define DEGREES (int) 5
#define TRUE (int) 1
#define FALSE (int) 0

typedef struct {
    float threshold;                                
    float default_value;                                
//    long warning_channel;                                
//    int  *warning_check_enable;                                
} DIV_STRUCT;

typedef struct {
    DIV_STRUCT divisor[1];//[19];                                
} DIV_TAB_STRUCT;

extern DIV_TAB_STRUCT div_tab;

int gauss (float *fx_ptr, float *fy_ptr, int degrees, int no_st, float *coeff_ptr);
void funk (float x, float *y_ptr);
int pos_def_matrix (float *, int , int);
float check_divisor(float, DIV_STRUCT *);

#endif