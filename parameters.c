#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
 
#include "parameters.h"
 
/* ------------------------------------------------------------ constants */
#define T_SUN_K        5772.0       /* K                                  */
#define M_BOL_SUN      4.74         /* absolute bolometric magnitude      */
#define L_SUN_W        3.828e26     /* W                                  */
#define LOGG_SUN       4.438        /* log10 of solar gravity, cgs        */
#define LIFETIME_SUN   1.0e10       /* years                              */
#define PC_TO_LY       3.26156
#define WIEN_NM_K      2.897771955e6   /* Wien constant in nm*K           */
 
/* mass-luminosity breakpoints (in solar units) */
#define ML_BREAK_1_MASS   0.43
#define ML_BREAK_2_MASS   2.0
#define ML_BREAK_3_MASS   55.0
 
/* ----------------------------------------------------------- tables */
typedef struct { double t; double v; } Anchor;
/* t = class_idx*10 + subclass, e.g. G2 -> 42, K0 -> 50, B5 -> 15 */
 
/* effective temperature, K */
static const Anchor T_V[] = {
    {5,41000},{10,31400},{12,20600},{13,17000},{15,15700},{18,12500},
    {20,9700},{22,8770},{25,8080},{30,7220},{32,6810},{35,6510},
    {40,5920},{42,5770},{45,5660},{50,5280},{52,5040},{55,4410},
    {60,3850},{62,3550},{65,3030} };
static const Anchor T_III[] = {
    {5,39000},{10,29000},{15,15000},{20,10000},{25,8000},{30,7100},
    {35,6500},{40,5900},{45,5100},{50,4800},{55,4000},{60,3850},{65,3400} };
static const Anchor T_I[] = {
    {5,38000},{10,26000},{15,14000},{20,9700},{25,8100},{30,7500},
    {35,6600},{40,5500},{45,4800},{50,4400},{55,3900},{60,3700},{65,3300} };
 
/* absolute visual magnitude M_V */
static const Anchor MV_V[] = {
    {5,-5.7},{10,-4.0},{15,-1.2},{20,0.6},{25,1.9},{30,2.7},{35,3.5},
    {40,4.4},{45,5.1},{50,5.9},{55,7.1},{60,8.8},{65,11.8} };
static const Anchor MV_III[] = {
    {5,-5.5},{10,-5.0},{15,-2.2},{20,-0.5},{25,0.5},{30,1.0},{35,0.9},
    {40,0.9},{45,0.7},{50,0.7},{55,-0.2},{60,-0.4},{65,-0.6} };
static const Anchor MV_I[] = {          /* crude mean of Ia and Ib */
    {5,-6.5},{10,-6.3},{15,-6.0},{20,-5.8},{25,-5.6},{30,-5.5},{35,-5.5},
    {40,-5.5},{45,-5.5},{50,-5.5},{55,-5.6},{60,-5.6},{65,-5.6} };
 
/* bolometric correction against temperature (ascending T) */
static const Anchor BC_T[] = {
    {3000,-2.6},{3500,-1.5},{4000,-0.9},{4500,-0.5},{5000,-0.25},
    {5500,-0.14},{6000,-0.07},{7000,0.0},{8000,-0.1},{9000,-0.2},
    {10000,-0.4},{15000,-1.4},{20000,-2.0},{30000,-3.0},{40000,-3.8} };
 
#define COUNT(a) ((int)(sizeof (a) / sizeof (a)[0]))
 
/* linear interpolation in an ascending table, clamped at both ends */
static double interp(const Anchor *a, int n, double x)
{
    if (x <= a[0].t)     return a[0].v;
    if (x >= a[n-1].t)   return a[n-1].v;
    for (int i = 1; i < n; i++) {
        if (x <= a[i].t) {
            double f = (x - a[i-1].t) / (a[i].t - a[i-1].t);
            return a[i-1].v + f * (a[i].v - a[i-1].v);
        }
    }
    return a[n-1].v;
}
 
/* value for the star's luminosity class; IV and II are the mean of neighbours */
static double by_lum_class(LumClass lum, double t,
                           const Anchor *v,   int nv,
                           const Anchor *iii, int niii,
                           const Anchor *i,   int ni)
{
    double dv = interp(v, nv, t), diii = interp(iii, niii, t), di = interp(i, ni, t);
    switch (lum) {
        case LUM_I:   return di;
        case LUM_II:  return 0.5 * (diii + di);
        case LUM_III: return diii;
        case LUM_IV:  return 0.5 * (dv + diii);
        default:      return dv;
    }
}
SpectralType parse_spectral(const char *s)
{
    SpectralType r;
    static const char CLASSES[] = "OBAFGKM";
    r.valid = 0; r.class_idx = -1; r.subclass = 0.0;
    r.lum = LUM_UNKNOWN; r.lum_assumed = 0; r.t = 0.0;
 
    if (s == NULL || s[0] == '\0') return r;
    const char *c = strchr(CLASSES, s[0]);
    if (c == NULL) return r;                 /* C, S, LBV, ... not supported */
 
    int idx = (int)(c - CLASSES);
    const char *p = s + 1;
 
    double sub = 0.0;
    while (isdigit((unsigned char)*p)) { sub = sub * 10.0 + (*p - '0'); p++; }
    if (*p == '.') {
        double scale = 0.1;
        p++;
        while (isdigit((unsigned char)*p)) { sub += (*p - '0') * scale; scale *= 0.1; p++; }
    }
 
    /* skip things like '+', '-', ':', '/' until the luminosity class */
    while (*p != '\0' && *p != 'I' && *p != 'V') p++;
 
    char run[4] = {0};
    int n = 0;
    while ((*p == 'I' || *p == 'V') && n < 3) run[n++] = *p++;
 
    if      (strcmp(run, "I")   == 0) r.lum = LUM_I;
    else if (strcmp(run, "II")  == 0) r.lum = LUM_II;
    else if (strcmp(run, "III") == 0) r.lum = LUM_III;
    else if (strcmp(run, "IV")  == 0) r.lum = LUM_IV;
    else if (strcmp(run, "V")   == 0) r.lum = LUM_V;
    else { r.lum = LUM_V; r.lum_assumed = 1; }   /* no class found */
 
    r.valid     = 1;
    r.class_idx = idx;
    r.subclass  = sub;
    r.t         = idx * 10.0 + sub;
    return r;
}
 
/* --------------------------------------------------------- small helpers */
static double temp_from_bv(double bv)
{
    if (!isfinite(bv) || bv < -0.4 || bv > 2.0) return NAN;
    return 4600.0 * (1.0 / (0.92 * bv + 1.7) + 1.0 / (0.92 * bv + 0.62));
}
 
/* invert the piecewise mass-luminosity relation (main sequence only) */
static double mass_from_lum(double lum)
{
    double l1 = pow(ML_BREAK_1_MASS, 4.0);
    double l2 = pow(ML_BREAK_2_MASS, 4.0);
    double l3 = 1.4 * pow(ML_BREAK_3_MASS, 3.5);
    if (lum <= 0.0)  return NAN;
    if (lum < l1)    return pow(lum / 0.23, 1.0 / 2.3);
    if (lum < l2)    return pow(lum, 1.0 / 4.0);
    if (lum < l3)    return pow(lum / 1.4, 1.0 / 3.5);
    return lum / 32000.0;
}
 
static void params_init(StarParams *p)
{
    p->spec.valid = 0; p->spec.class_idx = -1; p->spec.subclass = 0.0;
    p->spec.lum = LUM_UNKNOWN; p->spec.lum_assumed = 0; p->spec.t = 0.0;
    p->temp_K = p->dist_pc = p->dist_ly = p->dist_modulus = NAN;
    p->abs_mag = p->bol_corr = p->bol_mag = NAN;
    p->lum_solar = p->lum_watts = p->radius_solar = NAN;
    p->mass_solar = p->lifetime_yr = p->log_g = p->wien_nm = NAN;
    p->temp_method = p->dist_method = p->mass_method = METHOD_NONE;
}
 
/* ------------------------------------------------------- main function */
CalcStatus calculation_parameter(const Star *star, StarParams *out)
{
    if (star == NULL || out == NULL) return CALC_ERR_NULL_ARG;
    params_init(out);
    out->spec = parse_spectral(star->spectral);
 
    /* 1. temperature: stored value, then B-V, then spectral type */
    if (star->thermal_tag == THERMAL_TEMP_K &&
        isfinite(star->thermal.temp_K) && star->thermal.temp_K > 0.0) {
        out->temp_K = star->thermal.temp_K;
        out->temp_method = METHOD_GIVEN_TEMP;
    }
    if (!isfinite(out->temp_K) && star->thermal_tag == THERMAL_BV) {
        out->temp_K = temp_from_bv(star->thermal.bv);
        if (isfinite(out->temp_K)) out->temp_method = METHOD_BV;
    }
    if (!isfinite(out->temp_K) && out->spec.valid) {
        out->temp_K = by_lum_class(out->spec.lum, out->spec.t,
                                   T_V, COUNT(T_V), T_III, COUNT(T_III), T_I, COUNT(T_I));
        out->temp_method = METHOD_SPECTRAL;
    }
 
    /* 2. distance and absolute magnitude: parallax first, spectral type second */
    if (isfinite(star->parallax_mas) && star->parallax_mas > 0.0) {
        out->dist_pc  = 1000.0 / star->parallax_mas;
        out->abs_mag  = star->vmag - 5.0 * log10(out->dist_pc / 10.0);
        out->dist_method = METHOD_PARALLAX;
    } else if (out->spec.valid) {
        out->abs_mag  = by_lum_class(out->spec.lum, out->spec.t,
                                     MV_V, COUNT(MV_V), MV_III, COUNT(MV_III), MV_I, COUNT(MV_I));
        out->dist_pc  = 10.0 * pow(10.0, (star->vmag - out->abs_mag) / 5.0);
        out->dist_method = METHOD_SPECTRAL;
    }
    if (isfinite(out->dist_pc)) {
        out->dist_ly      = out->dist_pc * PC_TO_LY;
        out->dist_modulus = star->vmag - out->abs_mag;
    }
 
    /* 3. luminosity and radius (need both temperature and absolute magnitude) */
    if (isfinite(out->temp_K) && isfinite(out->abs_mag)) {
        out->bol_corr     = interp(BC_T, COUNT(BC_T), out->temp_K);
        out->bol_mag      = out->abs_mag + out->bol_corr;
        out->lum_solar    = pow(10.0, -0.4 * (out->bol_mag - M_BOL_SUN));
        out->lum_watts    = out->lum_solar * L_SUN_W;
        out->radius_solar = sqrt(out->lum_solar) * pow(T_SUN_K / out->temp_K, 2.0);
    }
    if (isfinite(out->temp_K)) out->wien_nm = WIEN_NM_K / out->temp_K;
 
    /* 4. mass: main-sequence stars only. Giants and supergiants stay NAN,
     *    because the mass-luminosity law does not apply to them. */
    if (isfinite(out->lum_solar) && out->spec.valid && out->spec.lum == LUM_V) {
        out->mass_solar  = mass_from_lum(out->lum_solar);
        out->mass_method = METHOD_MS_RELATION;
        if (isfinite(out->mass_solar)) {
            out->lifetime_yr = LIFETIME_SUN * out->mass_solar / out->lum_solar;
            if (isfinite(out->radius_solar))
                out->log_g = LOGG_SUN + log10(out->mass_solar)
                             - 2.0 * log10(out->radius_solar);
        }
    }
 
    if (!isfinite(out->temp_K) && !isfinite(out->dist_pc)) return CALC_ERR_NO_DATA;
    if (isfinite(out->lum_solar) && isfinite(out->radius_solar)) return CALC_OK;
    return CALC_PARTIAL;
}
 
/* ------------------------------------------------------------- report */
static const char *method_name(Method m)
{
    switch (m) {
        case METHOD_PARALLAX:    return "parallax";
        case METHOD_SPECTRAL:    return "spectral-type table";
        case METHOD_BV:          return "B-V colour";
        case METHOD_GIVEN_TEMP:  return "catalogue value";
        case METHOD_MS_RELATION: return "main-sequence relation";
        default:                 return "-";
    }
}
 
static void line(const char *label, double v, const char *unit, const char *how)
{
    if (isfinite(v)) printf("  %-24s %12.5g %-6s %s\n", label, v, unit, how);
    else             printf("  %-24s %12s\n", label, "n/a");
}
 
void print_star_report(const Star *star, const StarParams *p)
{
    const char *con = constellation_full(star->constellation);
    printf("\nHR %d  %s   constellation: %s\n", star->hr, star->name, con);
    printf("  spectral type (catalogue): %s   V = %.2f\n", star->spectral, star->vmag);
    printf("  RA %.4f h   Dec %+.4f deg (J2000)\n", star->ra_h, star->dec_deg);
    line("distance",          p->dist_pc,       "pc",  method_name(p->dist_method));
    line("distance",          p->dist_ly,       "ly",  "");
    line("absolute magnitude",p->abs_mag,       "mag", "");
    line("temperature",       p->temp_K,        "K",   method_name(p->temp_method));
    line("luminosity",        p->lum_solar,     "Lsun","");
    line("radius",            p->radius_solar,  "Rsun","");
    line("peak wavelength",   p->wien_nm,       "nm",  "");
    if (isfinite(p->mass_solar)) {
        line("mass",          p->mass_solar,    "Msun",method_name(p->mass_method));
        line("lifetime",      p->lifetime_yr,   "yr",  "");
        line("log g (cgs)",   p->log_g,         "",    "");
    } else {
        printf("  %-24s not estimated (not a main-sequence star or type unknown)\n", "mass");
    }
}
 


