#ifndef SKY_POSITION_H
#define SKY_POSITION_H
 
#include "data_access.h"
 
/* Everything the observer needs for one star at one instant and place. */
typedef struct {
    double lst_deg;        /* local sidereal time, 0..360          */
    double hour_angle_deg; /* 0..360, 0 = on the meridian          */
    double altitude_deg;   /* -90..+90, geometric (no refraction)  */
    double azimuth_deg;    /* 0..360 from North through East       */
    double airmass;        /* NAN when the star is below horizon   */
} SkyPosition;
 
/* ---- time ------------------------------------------------------------- */
/* UT calendar date -> Julian Date. Valid 1901..2099. Returns NAN if invalid. */
double julian_date(int year, int month, int day, double ut_hours);
 
/* Local sidereal time in DEGREES (0..360). Longitude east positive.
 * Divide by 15 for hours. */
double local_sidereal_time_deg(double jd, double longitude_east_deg);
 
/* ---- position --------------------------------------------------------- */
double normalize_deg(double deg);                     /* -> [0, 360)        */
 
/* hour angle = LST - RA, normalised to [0, 360). RA is given in hours. */
double hour_angle_deg(double lst_deg, double ra_h);
 
/* altitude and azimuth from latitude, declination and hour angle (degrees) */
void   altitude_azimuth(double lat_deg, double dec_deg, double ha_deg,
                        double *alt_deg, double *az_deg);
 
/* relative airmass for an altitude in degrees; NAN if altitude <= 0 */
double airmass_from_altitude(double alt_deg);
 
/* one call that does the whole chain for a catalogue star.
 * Returns 0 on success, -1 on NULL pointer, -2 on invalid date/latitude/longitude. */
int    star_sky_position(const Star *star, double jd,
                         double lat_deg, double lon_east_deg, SkyPosition *out);

#endif