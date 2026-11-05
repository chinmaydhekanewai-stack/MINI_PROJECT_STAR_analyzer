#ifndef PARAMETERS_H
#define PARAMETERS_H
 
#include "data_access.h"

/* ---------------------------------------------------------------- status */
typedef enum {
    CALC_OK,               /* all main parameters computed                   */
    CALC_PARTIAL,          /* some values are NAN (see the method fields)    */
    CALC_ERR_NULL_ARG,
    CALC_ERR_NO_DATA       /* neither temperature nor distance obtainable    */
} CalcStatus;
 
/* how a value was obtained, so the report can say so */
typedef enum {
    METHOD_NONE,
    METHOD_PARALLAX,       /* distance from catalogue parallax               */
    METHOD_SPECTRAL,       /* from spectral type tables                      */
    METHOD_BV,             /* temperature from B-V colour index              */
    METHOD_GIVEN_TEMP,     /* temperature stored directly in the catalogue   */
    METHOD_MS_RELATION     /* mass from the main-sequence mass-luminosity law*/
} Method;
 
typedef enum { LUM_UNKNOWN, LUM_V, LUM_IV, LUM_III, LUM_II, LUM_I } LumClass;
 
/* result of splitting a string such as "K0III" */
typedef struct {
    int      valid;        /* 0 for types the tables cannot handle (C, S, LBV) */
    int      class_idx;    /* 0=O 1=B 2=A 3=F 4=G 5=K 6=M                      */
    double   subclass;     /* 0.0 .. 9.5                                       */
    LumClass lum;
    int      lum_assumed;  /* 1 if no luminosity class was found (V assumed)   */
    double   t;            /* class_idx*10 + subclass : position in the tables */
} SpectralType;
 
/* every computed value for one star; unavailable values are NAN */
typedef struct {
    SpectralType spec;
 
    double temp_K;         Method temp_method;
    double dist_pc;        Method dist_method;
    double dist_ly;
    double dist_modulus;   /* m - M                                  */
    double abs_mag;        /* absolute visual magnitude M_V          */
    double bol_corr;       /* bolometric correction                  */
    double bol_mag;        /* M_bol                                  */
    double lum_solar;      /* L / L_sun                              */
    double lum_watts;      /* L in watts                             */
    double radius_solar;   /* R / R_sun                              */
    double mass_solar;     Method mass_method;
    double lifetime_yr;    /* main-sequence lifetime estimate        */
    double log_g;          /* log10 of surface gravity, cgs          */
    double wien_nm;        /* peak wavelength of black-body curve    */
} StarParams;
 
SpectralType parse_spectral(const char *s);
 
/* fills *out for one catalogue star; never modifies the catalogue */
CalcStatus   calculation_parameter(const Star *star, StarParams *out);
 
void         print_star_report(const Star *star, const StarParams *p);
 
#endif
 


