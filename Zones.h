#ifndef ZONES_H
#define ZONES_H


typedef enum {
    ZONE_OK,
    ZONE_ERR_NULL_ARG,
    ZONE_ERR_NOT_FOUND       /* city name or zone id not in the tables */
} ZoneStatus;

typedef struct {
    int         id;          /* 1..ZONE_COUNT, stable, used in request files */
    const char *name;
    double      lat_deg;     /* centre latitude  */
    double      lon_deg;     /* centre longitude */
} Zone;

typedef struct {
    const char *city;
    const char *state;
    double      lat_deg, lon_deg;   /* used only to pick the zone and to measure error */
} Place;

typedef struct {
    double dec_min_visible;  /* below this Dec a star never gets high enough  */
    double dec_max_visible;
    double dec_circumpolar;  /* at or above this Dec a star never sets        */
} DecWindow;

/* defined in zones.c */
extern const Zone  ZONES[];
extern const int   ZONE_COUNT;
extern const Place PLACES[];
extern const int   PLACE_COUNT;
double zone_distance_km(double lat1, double lon1, double lat2, double lon2);

int    zone_nearest(double lat, double lon);

int    zone_of_city(const char *city);
int    zone_index_by_id(int id);


ZoneStatus zone_center_of_city(const char *city, double *lat_deg,
                               double *lon_deg, int *zone_index);

/* same, starting from a zone id */
ZoneStatus zone_center_of_id(int id, double *lat_deg, double *lon_deg);

/* declination band for one zone.
 * min_alt_deg: lowest altitude you still call observable (for example 20).
 * margin_deg:  widens the band for users near the zone edge (for example 2). */
DecWindow zone_dec_window(const Zone *z, double min_alt_deg, double margin_deg);

#endif

