/*
 * zones.c  -  StellarCalc: India observing zones (see zones.h)
 *
 * Zone centres were chosen by hand and are approximate. Build with
 * -DZONES_DEMO to get a program that prints the worst city-to-centre
 * distance per zone, which is the real error of the zone idea:
 *   gcc -std=c99 -Wall -Wextra -DZONES_DEMO zones.c -o zones -lm
 */
#include <ctype.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "zones.h"

#define PI 3.14159265358979323846
#define EARTH_RADIUS_KM 6371.0088

/* ---------------------------------------------------------- zone table */
const Zone ZONES[] = {
    {  1, "Ladakh",                          34.5, 77.5 },
    {  2, "Jammu & Kashmir",                 33.5, 75.0 },
    {  3, "Punjab & Chandigarh",             31.0, 75.6 },
    {  4, "Himachal & Uttarakhand",          31.0, 78.0 },
    {  5, "Delhi, Haryana & West UP",        28.6, 77.3 },
    {  6, "Rajasthan West",                  27.0, 72.5 },
    {  7, "Rajasthan East & South",          25.5, 75.3 },
    {  8, "Uttar Pradesh Central & East",    26.3, 81.8 },
    {  9, "Bihar",                           25.6, 85.5 },
    { 10, "Jharkhand",                       23.6, 85.5 },
    { 11, "West Bengal",                     23.0, 88.0 },
    { 12, "Sikkim & North Bengal",           27.0, 88.5 },
    { 13, "Assam & Meghalaya",               26.0, 91.8 },
    { 14, "Arunachal Pradesh",               27.5, 94.2 },
    { 15, "North-East Hills",                24.5, 93.0 },
    { 16, "Gujarat North & Kutch",           23.2, 71.5 },
    { 17, "Gujarat South & Saurashtra",      21.8, 71.8 },
    { 18, "Madhya Pradesh West",             23.0, 76.3 },
    { 19, "Madhya Pradesh East",             24.5, 79.2 },
    { 20, "Chhattisgarh",                    21.5, 82.0 },
    { 21, "Maharashtra West",                18.6, 74.0 },
    { 22, "Maharashtra East",                20.7, 78.0 },
    { 23, "Goa & North Karnataka",           15.6, 75.0 },
    { 24, "Karnataka South",                 12.8, 76.3 },
    { 25, "Telangana",                       17.6, 79.0 },
    { 26, "Andhra Pradesh Coast",            17.0, 81.8 },
    { 27, "Andhra Pradesh South",            14.6, 79.0 },
    { 28, "Tamil Nadu North",                12.2, 79.3 },
    { 29, "Tamil Nadu South",                10.4, 78.0 },
    { 30, "Kerala",                          10.0, 76.4 },
    { 31, "Odisha",                          20.5, 84.5 },
    { 32, "Andaman & Nicobar",               11.7, 92.7 },
    { 33, "Lakshadweep",                     10.6, 72.6 },
};
const int ZONE_COUNT = (int)(sizeof ZONES / sizeof ZONES[0]);

/* --------------------------------------------------------- place table */
const Place PLACES[] = {
    {"Delhi","Delhi",28.61,77.21},          {"Meerut","Uttar Pradesh",28.98,77.71},
    {"Agra","Uttar Pradesh",27.18,78.01},   {"Lucknow","Uttar Pradesh",26.85,80.95},
    {"Kanpur","Uttar Pradesh",26.45,80.35}, {"Varanasi","Uttar Pradesh",25.32,82.97},
    {"Gorakhpur","Uttar Pradesh",26.76,83.37},
    {"Chandigarh","Chandigarh",30.73,76.78},{"Ludhiana","Punjab",30.90,75.86},
    {"Amritsar","Punjab",31.63,74.87},      {"Shimla","Himachal Pradesh",31.10,77.17},
    {"Dehradun","Uttarakhand",30.32,78.03}, {"Haridwar","Uttarakhand",29.95,78.16},
    {"Srinagar","Jammu & Kashmir",34.08,74.80},{"Jammu","Jammu & Kashmir",32.73,74.86},
    {"Leh","Ladakh",34.17,77.58},
    {"Jaipur","Rajasthan",26.91,75.79},     {"Jodhpur","Rajasthan",26.24,73.02},
    {"Bikaner","Rajasthan",28.02,73.31},    {"Udaipur","Rajasthan",24.59,73.71},
    {"Ahmedabad","Gujarat",23.02,72.57},    {"Bhuj","Gujarat",23.25,69.67},
    {"Rajkot","Gujarat",22.30,70.80},       {"Surat","Gujarat",21.17,72.83},
    {"Mumbai","Maharashtra",19.08,72.88},   {"Pune","Maharashtra",18.52,73.86},
    {"Nashik","Maharashtra",19.99,73.79},   {"Kolhapur","Maharashtra",16.70,74.24},
    {"Aurangabad","Maharashtra",19.88,75.34},{"Nagpur","Maharashtra",21.15,79.09},
    {"Panaji","Goa",15.50,73.83},           {"Hubballi","Karnataka",15.36,75.12},
    {"Bengaluru","Karnataka",12.97,77.59},  {"Mysuru","Karnataka",12.30,76.64},
    {"Mangaluru","Karnataka",12.91,74.86},
    {"Hyderabad","Telangana",17.39,78.49},  {"Warangal","Telangana",17.97,79.60},
    {"Vijayawada","Andhra Pradesh",16.51,80.65},{"Visakhapatnam","Andhra Pradesh",17.69,83.22},
    {"Tirupati","Andhra Pradesh",13.63,79.42},{"Nellore","Andhra Pradesh",14.44,79.99},
    {"Kurnool","Andhra Pradesh",15.83,78.04},
    {"Chennai","Tamil Nadu",13.08,80.27},   {"Puducherry","Puducherry",11.93,79.83},
    {"Tiruchirappalli","Tamil Nadu",10.80,78.69},{"Madurai","Tamil Nadu",9.92,78.12},
    {"Coimbatore","Tamil Nadu",11.02,76.96},
    {"Thiruvananthapuram","Kerala",8.52,76.94},{"Kochi","Kerala",9.93,76.27},
    {"Bhopal","Madhya Pradesh",23.26,77.41},{"Indore","Madhya Pradesh",22.72,75.86},
    {"Gwalior","Madhya Pradesh",26.22,78.18},{"Jabalpur","Madhya Pradesh",23.18,79.94},
    {"Raipur","Chhattisgarh",21.25,81.63},
    {"Bhubaneswar","Odisha",20.30,85.82},
    {"Patna","Bihar",25.59,85.14},          {"Gaya","Bihar",24.80,85.00},
    {"Ranchi","Jharkhand",23.34,85.31},     {"Jamshedpur","Jharkhand",22.80,86.20},
    {"Kolkata","West Bengal",22.57,88.36},  {"Siliguri","West Bengal",26.73,88.40},
    {"Gangtok","Sikkim",27.33,88.61},
    {"Guwahati","Assam",26.14,91.74},       {"Dibrugarh","Assam",27.48,94.91},
    {"Shillong","Meghalaya",25.57,91.88},   {"Itanagar","Arunachal Pradesh",27.08,93.61},
    {"Kohima","Nagaland",25.67,94.11},      {"Imphal","Manipur",24.82,93.95},
    {"Aizawl","Mizoram",23.73,92.72},       {"Agartala","Tripura",23.83,91.28},
    {"Port Blair","Andaman & Nicobar",11.62,92.73},
    {"Kavaratti","Lakshadweep",10.57,72.64},
};
const int PLACE_COUNT = (int)(sizeof PLACES / sizeof PLACES[0]);

/* ------------------------------------------------------------- helpers */
static double deg2rad(double d) { return d * PI / 180.0; }

static int equal_ignore_case(const char *a, const char *b)
{
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return 0;
        a++; b++;
    }
    return *a == *b;
}

/* -------------------------------------------------------- public API */
double zone_distance_km(double lat1, double lon1, double lat2, double lon2)
{
    double p1 = deg2rad(lat1), p2 = deg2rad(lat2);
    double dp = p2 - p1, dl = deg2rad(lon2 - lon1);
    double a = sin(dp / 2) * sin(dp / 2) +
               cos(p1) * cos(p2) * sin(dl / 2) * sin(dl / 2);
    if (a > 1.0) a = 1.0;
    return 2.0 * EARTH_RADIUS_KM * asin(sqrt(a));
}

int zone_nearest(double lat, double lon)
{
    int best = 0;
    double best_d = 1e30;
    for (int i = 0; i < ZONE_COUNT; i++) {
        double d = zone_distance_km(lat, lon, ZONES[i].lat_deg, ZONES[i].lon_deg);
        if (d < best_d) { best_d = d; best = i; }
    }
    return best;
}

int zone_of_city(const char *city)
{
    if (city == NULL) return -1;
    for (int i = 0; i < PLACE_COUNT; i++)
        if (equal_ignore_case(city, PLACES[i].city))
            return zone_nearest(PLACES[i].lat_deg, PLACES[i].lon_deg);
    return -1;
}

int zone_index_by_id(int id)
{
    for (int i = 0; i < ZONE_COUNT; i++)
        if (ZONES[i].id == id) return i;
    return -1;
}

ZoneStatus zone_center_of_city(const char *city, double *lat_deg,
                               double *lon_deg, int *zone_index)
{
    if (city == NULL || lat_deg == NULL || lon_deg == NULL) return ZONE_ERR_NULL_ARG;
    int z = zone_of_city(city);
    if (z < 0) return ZONE_ERR_NOT_FOUND;
    *lat_deg = ZONES[z].lat_deg;
    *lon_deg = ZONES[z].lon_deg;
    if (zone_index != NULL) *zone_index = z;
    return ZONE_OK;
}

ZoneStatus zone_center_of_id(int id, double *lat_deg, double *lon_deg)
{
    if (lat_deg == NULL || lon_deg == NULL) return ZONE_ERR_NULL_ARG;
    int z = zone_index_by_id(id);
    if (z < 0) return ZONE_ERR_NOT_FOUND;
    *lat_deg = ZONES[z].lat_deg;
    *lon_deg = ZONES[z].lon_deg;
    return ZONE_OK;
}

DecWindow zone_dec_window(const Zone *z, double min_alt_deg, double margin_deg)
{
    DecWindow w;
    double reach = 90.0 - min_alt_deg;          /* max |lat - dec| allowed */
    w.dec_min_visible = z->lat_deg - reach - margin_deg;
    w.dec_max_visible = z->lat_deg + reach + margin_deg;
    if (w.dec_min_visible < -90.0) w.dec_min_visible = -90.0;
    if (w.dec_max_visible >  90.0) w.dec_max_visible =  90.0;
    w.dec_circumpolar = 90.0 - z->lat_deg + margin_deg;
    return w;
}

/* ---------------------------------------------------------------- demo */
#ifdef ZONES_DEMO
int main(void)
{
    double worst[64] = {0};
    const char *worst_city[64] = {0};
    double overall = 0.0, sum = 0.0;

    for (int i = 0; i < PLACE_COUNT; i++) {
        int z = zone_nearest(PLACES[i].lat_deg, PLACES[i].lon_deg);
        double d = zone_distance_km(PLACES[i].lat_deg, PLACES[i].lon_deg,
                                    ZONES[z].lat_deg, ZONES[z].lon_deg);
        sum += d;
        if (d > overall) overall = d;
        if (d > worst[z]) { worst[z] = d; worst_city[z] = PLACES[i].city; }
    }
    printf("%d zones, %d places\n\n", ZONE_COUNT, PLACE_COUNT);
    for (int z = 0; z < ZONE_COUNT; z++)
        printf("%-3d %-30s %6.1f %6.1f  %s (%.0f km)\n", ZONES[z].id, ZONES[z].name,
               ZONES[z].lat_deg, ZONES[z].lon_deg,
               worst_city[z] ? worst_city[z] : "-", worst[z]);
    printf("\nMean city-to-centre distance: %.0f km\n", sum / PLACE_COUNT);
    printf("Worst city-to-centre distance: %.0f km (~%.1f deg of arc)\n",
           overall, overall / 111.2);
    return 0;
}
#endif