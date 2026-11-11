#include<gtest/gtest.h>

#include<cmath>
#include<cstring>
#include<string>

extern "C" {
#include "data_access.h"
#include "parameters.h"
#include "sky_position.h"
}

namespace {

const double kJ2000 = 2451545.0;

// A zero-initialised star with every optional field marked "unknown".
Star blank_star() {
    Star s;
    std::memset(&s, 0, sizeof s);
    s.name = "TEST";
    s.constellation = "";
    s.parallax_mas = NAN;
    s.thermal_tag = THERMAL_NONE;
    s.thermal.bv = NAN;
    return s;
}

} // End of anonymous namespace

// ============================== 1 ==========================================
TEST(Catalogue, FindStarByHr) {
    const Star* tau = find_star_by_hr(509);
    ASSERT_TRUE(tau != nullptr);
    EXPECT_EQ(tau->hr, 509);
    EXPECT_STREQ(tau->name, "Tau Ceti");

    EXPECT_TRUE(find_star_by_hr(-1) == nullptr);
    EXPECT_TRUE(find_star_by_hr(0) == nullptr);
    EXPECT_TRUE(find_star_by_hr(99999999) == nullptr);

    // every row must be found again, and must be the same row (pointer equality)
    ASSERT_GT(SEED_COUNT, 0);
    for (int i = 0; i < SEED_COUNT; i++) {
        SCOPED_TRACE(std::string("HR ") + std::to_string(SEED[i].hr));
        EXPECT_TRUE(find_star_by_hr(SEED[i].hr) == &SEED[i]);
    }
}

// ============================== 2 ==========================================
TEST(Catalogue, IntegrityOfAllRows) {
    ASSERT_GT(SEED_COUNT, 0);
    for (int i = 0; i < SEED_COUNT; i++) {
        const Star& s = SEED[i];
        SCOPED_TRACE(std::string("HR ") + std::to_string(s.hr));

        EXPECT_GT(s.hr, 0);
        ASSERT_TRUE(s.name != nullptr);
        EXPECT_GT(std::strlen(s.name), 0u);
        EXPECT_GE(s.ra_h, 0.0);
        EXPECT_LT(s.ra_h, 24.0);
        EXPECT_GE(s.dec_deg, -90.0);
        EXPECT_LE(s.dec_deg, 90.0);
        EXPECT_GT(s.vmag, -2.0);                    // Sirius is the brightest, -1.46
        EXPECT_LT(s.vmag, 8.0);                     // Bright Star Catalogue limit
        EXPECT_GT(std::strlen(s.spectral), 0u);
        EXPECT_TRUE(std::isnan(s.parallax_mas) || s.parallax_mas > 0.0);

        // a filled constellation code must exist in constellation_full()
        if (s.constellation != nullptr && s.constellation[0] != '\0') {
            EXPECT_STRNE(constellation_full(s.constellation), "unknown");
        }
        // HR numbers are primary keys: no duplicates
        for (int j = 0; j < i; j++) EXPECT_NE(SEED[j].hr, s.hr);
    }
    EXPECT_STREQ(constellation_full("Per"), "Perseus");
    EXPECT_STREQ(constellation_full("UMa"), "Ursa Major");
    EXPECT_STREQ(constellation_full("XXX"), "unknown");
    EXPECT_STREQ(constellation_full(""), "unknown");
}

// ============================== 3 ==========================================
TEST(Parameters, ParseSpectralTypes) {
    struct Case { const char* text; int valid; int cls; double sub; LumClass lum; int assumed; };
    const Case cases[] = {
        {"K0III",   1, 5, 0.0, LUM_III, 0},
        {"M0+II",   1, 6, 0.0, LUM_II,  0},
        {"B8IVp",   1, 1, 8.0, LUM_IV,  0},
        {"F5IV-V",  1, 3, 5.0, LUM_IV,  0},   // first class of a range wins
        {"A0Vvar",  1, 2, 0.0, LUM_V,   0},
        {"G9.5I",   1, 4, 9.5, LUM_I,   0},   // decimal subclass
        {"O9.7Ib",  1, 0, 9.7, LUM_I,   0},
        {"B0.5Iav", 1, 1, 0.5, LUM_I,   0},
        {"B3V+A",   1, 1, 3.0, LUM_V,   0},   // companion suffix ignored
        {"A0p",     1, 2, 0.0, LUM_V,   1},   // no class given: V assumed and flagged
    };
    for (const Case& c : cases) {
        SCOPED_TRACE(c.text);
        SpectralType t = parse_spectral(c.text);
        EXPECT_EQ(t.valid, c.valid);
        EXPECT_EQ(t.class_idx, c.cls);
        EXPECT_NEAR(t.subclass, c.sub, 1e-9);
        EXPECT_EQ(t.lum, c.lum);
        EXPECT_EQ(t.lum_assumed, c.assumed);
        EXPECT_NEAR(t.t, c.cls * 10.0 + c.sub, 1e-9);
    }
    // types the tables cannot handle must be flagged invalid, not guessed
    const char* invalid[] = {"LBV", "C6II", "S3.5/", "", "Wolf", "k0iii"};
    for (const char* s : invalid) {
        SCOPED_TRACE(s);
        EXPECT_EQ(parse_spectral(s).valid, 0);
    }
    EXPECT_EQ(parse_spectral(nullptr).valid, 0);
}

// ============================== 4 ==========================================
TEST(Parameters, TauCetiPhysicalSanity) {
    const Star* s = find_star_by_hr(509);          // G8V, a Sun-like dwarf
    ASSERT_TRUE(s != nullptr);
    StarParams p;
    EXPECT_EQ(calculation_parameter(s, &p), CALC_OK);

    EXPECT_EQ(p.temp_method, METHOD_SPECTRAL);
    EXPECT_EQ(p.dist_method, METHOD_SPECTRAL);
    EXPECT_GT(p.temp_K, 5000.0);   EXPECT_LT(p.temp_K, 5800.0);
    EXPECT_GT(p.dist_pc, 2.0);     EXPECT_LT(p.dist_pc, 6.0);    // real value ~3.65 pc
    EXPECT_GT(p.lum_solar, 0.3);   EXPECT_LT(p.lum_solar, 0.9);
    EXPECT_GT(p.radius_solar, 0.5); EXPECT_LT(p.radius_solar, 1.2);
    EXPECT_EQ(p.mass_method, METHOD_MS_RELATION);
    EXPECT_GT(p.mass_solar, 0.6);  EXPECT_LT(p.mass_solar, 1.1);
    EXPECT_GT(p.lifetime_yr, 5e9);
}

// ============================== 5 ==========================================
TEST(Parameters, SelfConsistencyAllStars) {
    int ok = 0;
    for (int i = 0; i < SEED_COUNT; i++) {
        const Star& s = SEED[i];
        SCOPED_TRACE(std::string("HR ") + std::to_string(s.hr));
        StarParams p;
        CalcStatus st = calculation_parameter(&s, &p);
        if (st != CALC_OK) continue;
        ok++;

        // Stefan-Boltzmann identity: L/Lsun = (R/Rsun)^2 * (T/Tsun)^4
        double rebuilt = p.radius_solar * p.radius_solar * std::pow(p.temp_K / 5772.0, 4.0);
        EXPECT_NEAR(rebuilt / p.lum_solar, 1.0, 1e-9);
        // unit and magnitude identities
        EXPECT_NEAR(p.dist_ly / p.dist_pc, 3.26156, 1e-9);
        EXPECT_NEAR(s.vmag - p.dist_modulus, p.abs_mag, 1e-9);
        EXPECT_NEAR(p.lum_watts / p.lum_solar, 3.828e26, 3.828e26 * 1e-9);
        // physical ranges
        EXPECT_GT(p.temp_K, 1500.0);   EXPECT_LT(p.temp_K, 60000.0);
        EXPECT_GT(p.dist_pc, 0.0);
        EXPECT_GT(p.lum_solar, 0.0);
        EXPECT_GT(p.radius_solar, 0.0);
        EXPECT_GT(p.wien_nm, 0.0);
    }
    EXPECT_GE(ok, SEED_COUNT / 2);   // most stars must be computable
}

// ============================== 6 ==========================================
TEST(Parameters, StatusCodesAndMassRule) {
    StarParams p;
    const Star* tau = find_star_by_hr(509);
    ASSERT_TRUE(tau != nullptr);
    EXPECT_EQ(calculation_parameter(nullptr, &p), CALC_ERR_NULL_ARG);
    EXPECT_EQ(calculation_parameter(tau, nullptr), CALC_ERR_NULL_ARG);

    // types the tables cannot handle, and no parallax/B-V: nothing can be computed
    const int no_data[] = {4210 /* Eta Carinae, LBV */, 423 /* carbon star */};
    for (int hr : no_data) {
        SCOPED_TRACE(std::string("HR ") + std::to_string(hr));
        const Star* s = find_star_by_hr(hr);
        ASSERT_TRUE(s != nullptr);
        EXPECT_EQ(calculation_parameter(s, &p), CALC_ERR_NO_DATA);
        EXPECT_TRUE(std::isnan(p.temp_K));
        EXPECT_TRUE(std::isnan(p.lum_solar));
    }

    // the main-sequence mass law must NOT be applied to a supergiant
    const Star* betel = find_star_by_hr(2061);
    ASSERT_TRUE(betel != nullptr);
    EXPECT_EQ(calculation_parameter(betel, &p), CALC_OK);
    EXPECT_TRUE(std::isnan(p.mass_solar));
    EXPECT_EQ(p.mass_method, METHOD_NONE);
    EXPECT_GT(p.radius_solar, 100.0);               // a very large star

    // a hotter star of the same size class is more luminous: Sirius A1V vs Tau Ceti G8V
    StarParams sirius, tauc;
    ASSERT_TRUE(find_star_by_hr(2491) != nullptr);
    calculation_parameter(find_star_by_hr(2491), &sirius);
    calculation_parameter(tau, &tauc);
    EXPECT_GT(sirius.temp_K, tauc.temp_K);
    EXPECT_GT(sirius.lum_solar, tauc.lum_solar);
}

// ============================== 7 ==========================================
TEST(Parameters, ParallaxAndBvRoutes) {
    // Sun-like star at exactly 10 pc (parallax 100 mas), temperature given directly
    Star s = blank_star();
    s.vmag = 4.83;
    s.parallax_mas = 100.0;
    std::strcpy(s.spectral, "G2V");
    s.thermal_tag = THERMAL_TEMP_K;
    s.thermal.temp_K = 5772.0;

    StarParams p;
    EXPECT_EQ(calculation_parameter(&s, &p), CALC_OK);
    EXPECT_EQ(p.dist_method, METHOD_PARALLAX);
    EXPECT_EQ(p.temp_method, METHOD_GIVEN_TEMP);
    EXPECT_NEAR(p.dist_pc, 10.0, 1e-9);
    EXPECT_NEAR(p.abs_mag, 4.83, 1e-9);             // at 10 pc, M = m
    EXPECT_NEAR(p.lum_solar, 1.0, 0.07);            // the Sun
    EXPECT_NEAR(p.radius_solar, 1.0, 0.04);
    EXPECT_NEAR(p.mass_solar, 1.0, 0.07);

    // B-V route: 0.65 is the Sun's colour, about 5800 K
    s.thermal_tag = THERMAL_BV;
    s.thermal.bv = 0.65;
    EXPECT_EQ(calculation_parameter(&s, &p), CALC_OK);
    EXPECT_EQ(p.temp_method, METHOD_BV);
    EXPECT_NEAR(p.temp_K, 5800.0, 150.0);

    // an impossible B-V falls back to the spectral type
    s.thermal.bv = 5.0;
    calculation_parameter(&s, &p);
    EXPECT_EQ(p.temp_method, METHOD_SPECTRAL);

    // zero or negative parallax must not be used for distance
    s.parallax_mas = -5.0;
    calculation_parameter(&s, &p);
    EXPECT_EQ(p.dist_method, METHOD_SPECTRAL);
    s.parallax_mas = 0.0;
    calculation_parameter(&s, &p);
    EXPECT_EQ(p.dist_method, METHOD_SPECTRAL);
}

// ============================== 8 ==========================================
TEST(SkyTime, JulianAndSiderealTime) {
    // Julian Date: known values
    EXPECT_NEAR(julian_date(2000, 1, 1, 12.0), kJ2000, 1e-9);
    EXPECT_NEAR(julian_date(2000, 1, 1, 0.0), 2451544.5, 1e-9);
    EXPECT_NEAR(julian_date(2020, 1, 1, 0.0), 2458849.5, 1e-9);
    EXPECT_NEAR(julian_date(2024, 3, 1, 0.0) - julian_date(2024, 2, 28, 0.0), 2.0, 1e-9);  // leap year
    EXPECT_NEAR(julian_date(2023, 3, 1, 0.0) - julian_date(2023, 2, 28, 0.0), 1.0, 1e-9);

    // invalid dates are rejected with NaN
    EXPECT_TRUE(std::isnan(julian_date(2023, 2, 29, 0.0)));   // not a leap year
    EXPECT_TRUE(std::isnan(julian_date(2026, 2, 30, 0.0)));
    EXPECT_TRUE(std::isnan(julian_date(2026, 13, 1, 0.0)));
    EXPECT_TRUE(std::isnan(julian_date(2026, 0, 1, 0.0)));
    EXPECT_TRUE(std::isnan(julian_date(2026, 4, 31, 0.0)));
    EXPECT_TRUE(std::isnan(julian_date(2026, 4, 0, 0.0)));
    EXPECT_TRUE(std::isnan(julian_date(1900, 1, 1, 0.0)));    // outside 1901..2099
    EXPECT_TRUE(std::isnan(julian_date(2100, 1, 1, 0.0)));
    EXPECT_TRUE(std::isnan(julian_date(2026, 1, 1, 24.0)));
    EXPECT_TRUE(std::isnan(julian_date(2026, 1, 1, -1.0)));
    EXPECT_FALSE(std::isnan(julian_date(2000, 2, 29, 0.0)));  // 2000 is a leap year

    // Local sidereal time: known values
    EXPECT_NEAR(local_sidereal_time_deg(kJ2000, 0.0), 280.46062, 1e-4);
    EXPECT_NEAR(local_sidereal_time_deg(2451544.5, 0.0), 99.9678, 1e-3);   // 6h39m52s at 2000-01-01 0h UT
    EXPECT_NEAR(local_sidereal_time_deg(kJ2000, 15.0), 295.46062, 1e-4);   // +15 deg east = +1 hour
    EXPECT_NEAR(local_sidereal_time_deg(kJ2000, 100.0), 20.46062, 1e-4);   // wraps past 360

    // output is always in [0, 360), for any date and longitude
    for (int k = 0; k < 60; k++) {
        for (int lon = -180; lon <= 180; lon += 45) {
            double lst = local_sidereal_time_deg(kJ2000 + k * 123.456 - 3000.0, lon);
            EXPECT_GE(lst, 0.0);
            EXPECT_LT(lst, 360.0);
        }
    }
    // one sidereal day later gives the same LST (to within the rate formula)
    double a = local_sidereal_time_deg(kJ2000 + 5000.25, 77.0);
    double b = local_sidereal_time_deg(kJ2000 + 5000.25 + 0.99726957, 77.0);
    double diff = std::fmod(std::fabs(a - b), 360.0);
    EXPECT_TRUE(diff < 0.01 || diff > 359.99);
}

// ============================== 9 ==========================================
TEST(SkyGeometry, KnownCasesAndInvariants) {
    double alt = 0.0, az = 0.0;

    // declination equal to latitude, on the meridian: passes through the zenith
    altitude_azimuth(16.7, 16.7, 0.0, &alt, &az);
    EXPECT_NEAR(alt, 90.0, 1e-6);

    // transit altitude = 90 - |lat - dec|; azimuth is South (180) or North (0)
    altitude_azimuth(28.6, 10.0, 0.0, &alt, &az);
    EXPECT_NEAR(alt, 71.4, 1e-6);
    EXPECT_NEAR(az, 180.0, 1e-6);
    altitude_azimuth(28.6, 40.0, 0.0, &alt, &az);
    EXPECT_NEAR(alt, 78.6, 1e-6);
    EXPECT_TRUE(az < 1e-4 || az > 360.0 - 1e-4);

    // before transit (east) azimuth < 180; after transit (west) azimuth > 180
    altitude_azimuth(16.7, -10.0, 330.0, &alt, &az);
    EXPECT_GT(az, 0.0);   EXPECT_LT(az, 180.0);
    altitude_azimuth(16.7, -10.0, 30.0, &alt, &az);
    EXPECT_GT(az, 180.0); EXPECT_LT(az, 360.0);

    // the celestial pole stands at an altitude equal to the latitude, at any hour angle
    for (int ha = 0; ha < 360; ha += 60) {
        altitude_azimuth(16.7, 90.0, ha, &alt, &az);
        EXPECT_NEAR(alt, 16.7, 1e-6);
    }
    // a star near the pole (Polaris, dec 89.26) stays within 0.8 deg of the latitude
    for (int ha = 0; ha < 360; ha += 30) {
        altitude_azimuth(28.6, 89.26, ha, &alt, &az);
        EXPECT_NEAR(alt, 28.6, 0.8);
    }
    // north/south mirror symmetry at the meridian
    double alt_n, alt_s;
    altitude_azimuth(20.0, 5.0, 0.0, &alt_n, &az);
    altitude_azimuth(-20.0, -5.0, 0.0, &alt_s, &az);
    EXPECT_NEAR(alt_n, alt_s, 1e-9);

    // hour angle: normalised to [0, 360)
    EXPECT_NEAR(hour_angle_deg(10.0, 1.0), 355.0, 1e-9);    // LST 10 deg, RA 1 h = 15 deg
    EXPECT_NEAR(hour_angle_deg(100.0, 6.0), 10.0, 1e-9);
    EXPECT_NEAR(hour_angle_deg(0.0, 0.0), 0.0, 1e-9);
    EXPECT_NEAR(normalize_deg(-30.0), 330.0, 1e-9);
    EXPECT_NEAR(normalize_deg(725.0), 5.0, 1e-9);

    // ranges hold across a full sweep of latitude, declination and hour angle
    for (int lat = -90; lat <= 90; lat += 30)
        for (int dec = -90; dec <= 90; dec += 30)
            for (int ha = 0; ha < 360; ha += 30) {
                altitude_azimuth(lat, dec, ha, &alt, &az);
                EXPECT_GE(alt, -90.0); EXPECT_LE(alt, 90.0);
                EXPECT_GE(az, 0.0);    EXPECT_LT(az, 360.0);
                EXPECT_FALSE(std::isnan(alt));
                EXPECT_FALSE(std::isnan(az));
            }
}

// ============================== 10 =========================================
TEST(SkyGeometry, AirmassAndFullChain) {
    // airmass: 1 at the zenith, grows toward the horizon, undefined below it
    EXPECT_NEAR(airmass_from_altitude(90.0), 1.0, 1e-3);
    EXPECT_NEAR(airmass_from_altitude(30.0), 2.0, 0.02);
    EXPECT_TRUE(std::isnan(airmass_from_altitude(0.0)));
    EXPECT_TRUE(std::isnan(airmass_from_altitude(-10.0)));
    for (int a = 5; a < 90; a += 5)
        EXPECT_GT(airmass_from_altitude(a), airmass_from_altitude(a + 5));

    // star_sky_position: input rejection
    Star s = blank_star();
    s.ra_h = 6.7525;
    s.dec_deg = -16.7161;
    SkyPosition pos;
    const double jd = julian_date(2026, 11, 11, 14.5);
    EXPECT_EQ(star_sky_position(nullptr, jd, 16.7, 74.2, &pos), -1);
    EXPECT_EQ(star_sky_position(&s, jd, 16.7, 74.2, nullptr), -1);
    EXPECT_EQ(star_sky_position(&s, jd, 95.0, 74.2, &pos), -2);
    EXPECT_EQ(star_sky_position(&s, jd, -91.0, 74.2, &pos), -2);
    EXPECT_EQ(star_sky_position(&s, jd, 16.7, 200.0, &pos), -2);
    EXPECT_EQ(star_sky_position(&s, NAN, 16.7, 74.2, &pos), -2);

    // valid call: fields in range, consistent with the step-by-step functions
    ASSERT_EQ(star_sky_position(&s, jd, 16.7, 74.2, &pos), 0);
    EXPECT_GE(pos.hour_angle_deg, 0.0);  EXPECT_LT(pos.hour_angle_deg, 360.0);
    EXPECT_GE(pos.azimuth_deg, 0.0);     EXPECT_LT(pos.azimuth_deg, 360.0);
    double alt, az;
    altitude_azimuth(16.7, s.dec_deg, pos.hour_angle_deg, &alt, &az);
    EXPECT_NEAR(pos.altitude_deg, alt, 1e-9);
    EXPECT_NEAR(pos.lst_deg, local_sidereal_time_deg(jd, 74.2), 1e-9);
    // a below-horizon star has no airmass
    if (pos.altitude_deg <= 0.0) EXPECT_TRUE(std::isnan(pos.airmass));

    // Polaris from the catalogue through the full chain: altitude stays within
    // 1 degree of the observer's latitude at every hour of the night
    const Star* polaris = find_star_by_hr(424);
    ASSERT_TRUE(polaris != nullptr);
    const double lat = 28.6, lon = 77.2;
    for (int h = 0; h < 24; h++) {
        SCOPED_TRACE(std::string("hour ") + std::to_string(h));
        SkyPosition p;
        ASSERT_EQ(star_sky_position(polaris, julian_date(2026, 11, 11, 0.0) + h / 24.0, lat, lon, &p), 0);
        EXPECT_NEAR(p.altitude_deg, lat, 1.0);
        EXPECT_FALSE(std::isnan(p.airmass));
    }
}