#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "data_access.h"
#include "parameters.h"
#include "sky_position.h"
#include "zones.h"
 
#define MIN_ALT_DEG    20.0    
#define IST_OFFSET_MIN 330.0   /* India Standard Time = UT + 5 h 30 min       */
#define MAX_LIST       512
 
/* ------------------------------------------------------------- input helpers */
static int read_line(const char *prompt, char *buf, size_t size)
{
    printf("%s", prompt);
    if (fgets(buf, (int)size, stdin) == NULL) return 0;
    buf[strcspn(buf, "\r\n")] = '\0';
    return 1;
}
 
/* whole-number parse: returns 1 on success */
static int parse_long(const char *text, long *out)
{
    char *end;
    if (text[0] == '\0') return 0;
    *out = strtol(text, &end, 10);
    return *end == '\0';
}

static const Star *ask_star(void)
{
    char line[64];
    long hr;
    if (!read_line("Enter HR number: ", line, sizeof line)) return NULL;
    if (!parse_long(line, &hr)) {
        printf("Invalid input: please enter a whole number.\n");
        return NULL;
    }
    const Star *s = find_star_by_hr((int)hr);
    if (s == NULL) printf("HR %ld is not in the catalogue.\n", hr);
    return s;
}
 
static void list_cities(void)
{
    printf("Available cities:\n");
    for (int i = 0; i < PLACE_COUNT; i++)
        printf("  %-20s%s", PLACES[i].city, (i % 3 == 2) ? "\n" : "");
    printf("\n");
}
 
/* asks for a city; returns 1 and fills lat/lon/zone on success, 0 if the user gives up */
static int ask_location(double *lat, double *lon, int *zone_index)
{
    char city[64];
    for (;;) {
        if (!read_line("City (or 'list', empty to cancel): ", city, sizeof city)) return 0;
        if (city[0] == '\0') return 0;
        if (strcmp(city, "list") == 0) { list_cities(); continue; }
        if (zone_center_of_city(city, lat, lon, zone_index) == ZONE_OK) return 1;
        printf("City '%s' not found. Type 'list' to see the cities.\n", city);
    }
}
static double ask_time_ist(void)
{
    char line[64];
    int y, mo, d, h, mi, used = 0;
 
    if (!read_line("Date and time in IST (YYYY-MM-DD HH:MM): ", line, sizeof line)) return NAN;
    if (sscanf(line, "%d-%d-%d %d:%d%n", &y, &mo, &d, &h, &mi, &used) != 5 ||
        line[used] != '\0') {
        printf("Invalid format. Example: 2026-11-11 20:00\n");
        return NAN;
    }
    if (h < 0 || h > 23 || mi < 0 || mi > 59) {
        printf("Invalid time: hour must be 0-23 and minute 0-59.\n");
        return NAN;
    }
    double jd0 = julian_date(y, mo, d, 0.0);
    if (isnan(jd0)) {
        printf("Invalid date (check day and month; supported years are 1901-2099).\n");
        return NAN;
    }
    return jd0 + ((h * 60.0 + mi) - IST_OFFSET_MIN) / 1440.0;
}
 
static const char *altitude_comment(double alt)
{
    if (alt <= 0.0)        return "below the horizon";
    if (alt < MIN_ALT_DEG) return "up, but low (poor viewing)";
    return "well placed for observing";
}
 
static const char *compass(double az)
{
    static const char *names[] = {"N", "NE", "E", "SE", "S", "SW", "W", "NW"};
    return names[(int)((az + 22.5) / 45.0) % 8];
}
 
/* --------------------------------------------------------------- menu actions */
static void action_star_info(void)
{
    const Star *star = ask_star();
    if (star == NULL) return;
 
    StarParams p;
    CalcStatus st = calculation_parameter(star, &p);
    if (st == CALC_ERR_NO_DATA || st == CALC_ERR_NULL_ARG) {
        printf("Cannot compute parameters for HR %d (status %d).\n", star->hr, (int)st);
        return;
    }
    print_star_report(star, &p);
    if (st == CALC_PARTIAL)
        printf("\n(note: some values are unavailable for this star)\n");
}
 
static void action_star_position(void)
{
    const Star *star = ask_star();
    if (star == NULL) return;
 
    double lat, lon;
    int zi;
    if (!ask_location(&lat, &lon, &zi)) return;
    double jd = ask_time_ist();
    if (isnan(jd)) return;
 
    SkyPosition pos;
    if (star_sky_position(star, jd, lat, lon, &pos) != 0) {
        printf("Could not compute the position (invalid input).\n");
        return;
    }
    printf("\nHR %d  %s\n", star->hr, star->name);
    printf("  zone %d: %s (centre %.1f N, %.1f E)\n",
           ZONES[zi].id, ZONES[zi].name, lat, lon);
    printf("  local sidereal time : %.2f h\n", pos.lst_deg / 15.0);
    printf("  altitude            : %+.1f deg\n", pos.altitude_deg);
    printf("  azimuth             : %.1f deg (%s)\n", pos.azimuth_deg, compass(pos.azimuth_deg));
    if (pos.altitude_deg > 0.0) printf("  airmass             : %.2f\n", pos.airmass);
    printf("  status              : %s\n", altitude_comment(pos.altitude_deg));
    printf("  accuracy            : about 1-2 deg (zone centre used, no refraction)\n");
}
 
typedef struct { int index; double alt; double az; } Row;
 
static int by_altitude_desc(const void *a, const void *b)
{
    double da = ((const Row *)a)->alt, db = ((const Row *)b)->alt;
    return (da < db) - (da > db);
}
 
static void action_whats_up(void)
{
    double lat, lon;
    int zi;
    if (!ask_location(&lat, &lon, &zi)) return;
    double jd = ask_time_ist();
    if (isnan(jd)) return;
 
    static Row rows[MAX_LIST];
    int n = 0;
    for (int i = 0; i < SEED_COUNT && n < MAX_LIST; i++) {
        SkyPosition pos;
        if (star_sky_position(&SEED[i], jd, lat, lon, &pos) != 0) continue;
        if (pos.altitude_deg >= MIN_ALT_DEG) {
            rows[n].index = i;
            rows[n].alt   = pos.altitude_deg;
            rows[n].az    = pos.azimuth_deg;
            n++;
        }
    }
    qsort(rows, (size_t)n, sizeof rows[0], by_altitude_desc);

    printf("\nZone %d: %s. Stars at least %.0f deg above the horizon:\n",
           ZONES[zi].id, ZONES[zi].name, MIN_ALT_DEG);
    if (n == 0) printf("  none in the catalogue at this time.\n");
    for (int k = 0; k < n; k++) {
        const Star *s = &SEED[rows[k].index];
        printf("  HR %-5d %-18s V %5.2f   alt %5.1f   az %5.1f %s\n",
               s->hr, s->name, s->vmag, rows[k].alt, rows[k].az, compass(rows[k].az));
    }
    printf("\n(note: the Sun is not modelled yet, so this list does not check for daylight)\n");
}
int main(void)
{
    char line[64];
 
    printf("StellarCalc - %d stars, %d zones (catalogue version %d)\n",
           SEED_COUNT, ZONE_COUNT, CATALOG_VERSION);
 
    for (;;) {
        long choice;
        printf("\n1) Star information\n2) Star position at a place and time\n"
               "3) Which stars are up at a place and time\n0) Exit\n");
        if (!read_line("Choice: ", line, sizeof line)) break;
        if (!parse_long(line, &choice)) {
            printf("Invalid input: please enter a number from the menu.\n");
            continue;
        }
        if (choice == 0) break;
        switch (choice) {
            case 1: action_star_info();     break;
            case 2: action_star_position(); break;
            case 3: action_whats_up();      break;
            default: printf("Please choose 0, 1, 2 or 3.\n");
        }
    }
    return 0;
}
 
