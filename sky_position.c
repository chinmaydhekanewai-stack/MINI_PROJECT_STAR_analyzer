#include <math.h>
#include <stddef.h>
 
#include "sky_position.h"
 
#define PI 3.14159265358979323846
#define DEG2RAD(d) ((d) * PI / 180.0)
#define RAD2DEG(r) ((r) * 180.0 / PI)
 
static double clamp(double x, double lo, double hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}
 
static int is_leap(int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; }
 
/* ------------------------------------------------------------------ time */
double julian_date(int year, int month, int day, double ut_hours)
{
    static const int dim[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    if (year < 1901 || year > 2099 || month < 1 || month > 12) return NAN;
    int maxd = dim[month - 1] + ((month == 2 && is_leap(year)) ? 1 : 0);
    if (day < 1 || day > maxd || ut_hours < 0.0 || ut_hours >= 24.0) return NAN;
 
    return 367.0 * year
         - floor(7.0 * (year + floor((month + 9) / 12.0)) / 4.0)
         + floor(275.0 * month / 9.0)
         + day + 1721013.5 + ut_hours / 24.0;
}
 
double normalize_deg(double deg)
{
    double r = fmod(deg, 360.0);
    return r < 0.0 ? r + 360.0 : r;    
}
 
double local_sidereal_time_deg(double jd, double longitude_east_deg)
{
    double n    = jd - 2451545.0;                       /* days since J2000.0 */
    double gmst = 280.46061837 + 360.98564736629 * n;   /* degrees            */
    return normalize_deg(gmst + longitude_east_deg);
}
 
/* -------------------------------------------------------------- position */
double hour_angle_deg(double lst_deg, double ra_h)
{
    return normalize_deg(lst_deg - 15.0 * ra_h);
}
 
void altitude_azimuth(double lat_deg, double dec_deg, double ha_deg,
                      double *alt_deg, double *az_deg)
{
    double lat = DEG2RAD(lat_deg), dec = DEG2RAD(dec_deg), ha = DEG2RAD(ha_deg);
 
    double sin_alt = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(ha);
    sin_alt = clamp(sin_alt, -1.0, 1.0);
    double alt = asin(sin_alt);
 
    double denom = cos(alt) * cos(lat);
    double az;
    if (fabs(denom) < 1e-12) {
        az = 0.0;                         /* zenith or a pole: azimuth undefined */
    } else {
        double cos_az = clamp((sin(dec) - sin_alt * sin(lat)) / denom, -1.0, 1.0);
        az = acos(cos_az);
        if (sin(ha) > 0.0) az = 2.0 * PI - az;   /* star is west of meridian */
    }
    *alt_deg = RAD2DEG(alt);
    *az_deg  = normalize_deg(RAD2DEG(az));
}
 
double airmass_from_altitude(double alt_deg)
{
    if (alt_deg <= 0.0) return NAN;
    return 1.0 / (sin(DEG2RAD(alt_deg)) + 0.50572 * pow(alt_deg + 6.07995, -1.6364));
}
 
int star_sky_position(const Star *star, double jd,
                      double lat_deg, double lon_east_deg, SkyPosition *out)
{
    if (star == NULL || out == NULL) return -1;
    if (!isfinite(jd) || lat_deg < -90.0 || lat_deg > 90.0 ||
        lon_east_deg < -180.0 || lon_east_deg > 180.0) return -2;
 
    out->lst_deg        = local_sidereal_time_deg(jd, lon_east_deg);
    out->hour_angle_deg = hour_angle_deg(out->lst_deg, star->ra_h);
    altitude_azimuth(lat_deg, star->dec_deg, out->hour_angle_deg,
                     &out->altitude_deg, &out->azimuth_deg);
    out->airmass        = airmass_from_altitude(out->altitude_deg);
    return 0;
}
 