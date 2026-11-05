#ifndef DATA_ACCESS_H
#define DATA_ACCESS_H
 
#define CATALOG_VERSION 2
 
typedef enum { REGION_MAIN_SEQUENCE, REGION_SUBGIANT, REGION_GIANT,
               REGION_SUPERGIANT, REGION_AGB_CARBON_S } HRRegion;
 
typedef enum { THERMAL_NONE, THERMAL_BV, THERMAL_TEMP_K } ThermalTag;
 
typedef struct {
    int         hr;              /* Harvard Revised number, primary key   */
    const char *name;            /* common name or "HR nnn"               */
    const char *constellation;   /* IAU 3-letter code, "" if unknown      */
    double      ra_h;            /* J2000, decimal hours                  */
    double      dec_deg;         /* J2000, decimal degrees                */
    double      pm_ra, pm_dec;   /* proper motion, catalogue units        */
    double      vmag;            /* apparent visual magnitude             */
    double      parallax_mas;    /* NAN until filled                      */
    char        spectral[12];    /* verbatim from catalogue               */
    HRRegion    region;          /* provisional                           */
    ThermalTag  thermal_tag;     /* says which union member is valid      */
    union { double bv; double temp_K; } thermal;
} Star;
 
/* defined in data_access.c */
extern const Star SEED[];
extern const int  SEED_COUNT;
 
/* IAU abbreviation -> full constellation name ("unknown" if not found) */
const char *constellation_full(const char *abbr);
 
/* returns a pointer into SEED, or NULL if no star has this HR number */
const Star *find_star_by_hr(int hr);
 
#endif
 