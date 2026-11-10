#include <stdio.h>
#include <stdlib.h>
 
#include "parameters.h"
#include "data_access.h"
 
int main(void)
{
    char line[64];
    long hr;
    char *end;

    printf("StellarCalc - %d stars in catalogue (version %d)\n", SEED_COUNT, CATALOG_VERSION);
    printf("Enter HR number: ");
 
    if (fgets(line, sizeof line, stdin) == NULL) {
        printf("No input.\n");
        return 1;
    }
    hr = strtol(line, &end, 10);
    if (end == line || (*end != '\n' && *end != '\0')) {
        printf("Invalid input: please enter a whole number.\n");
        return 1;
    }
 
    const Star *star = find_star_by_hr((int)hr);
    if (star == NULL) {
        printf("HR %ld is not in the catalogue.\n", hr);
        return 2;
    }
 
    StarParams p;
    CalcStatus st = calculation_parameter(star, &p);
    if (st == CALC_ERR_NO_DATA || st == CALC_ERR_NULL_ARG) {
        printf("Cannot compute parameters for HR %ld (status %d).\n", hr, (int)st);
        return 3;
    }
 
    print_star_report(star, &p);
    if (st == CALC_PARTIAL)
        printf("\n(note: some values are unavailable for this star)\n");
    return 0;
}
 