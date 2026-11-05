#include <math.h>
#include <stdio.h>

#define CATALOG_VERSION 2

typedef enum { REGION_MAIN_SEQUENCE, REGION_SUBGIANT, REGION_GIANT,
               REGION_SUPERGIANT, REGION_AGB_CARBON_S } HRRegion;

typedef enum { THERMAL_NONE, THERMAL_BV, THERMAL_TEMP_K } ThermalTag;

typedef struct {
    int         hr;              /* Harvard Revised number, primary key   */
    const char *name;            /* common name or "HR nnn"               */
    const char *constellation;   /* fill from catalogue                   */
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
    #define RA(h,m,s)        ((h) + (m)/60.0 + (s)/3600.0)
#define DEC(sg,d,m,s)    ((sg) * ((d) + (m)/60.0 + (s)/3600.0))

#define S(hr_,nm_,con_,ra_,dec_,pmr_,pmd_,v_,sp_,reg_) \
    { .hr=(hr_), .name=(nm_), .constellation=(con_), .ra_h=(ra_), .dec_deg=(dec_), \
    .pm_ra=(pmr_), .pm_dec=(pmd_), .vmag=(v_), .parallax_mas=NAN, \
    .spectral=(sp_), .region=(reg_), .thermal_tag=THERMAL_NONE, .thermal={.bv=NAN} }

static const Star SEED[] = {
 /* --- Original bright, well-known (V < 3) --- */
S(15,"Alpheratz","And",   RA(0, 8,23.3), DEC(+1,29, 5,26),  0.136,-0.163, 2.06,"B8IVp",  REGION_SUBGIANT),
S(21,"Caph","Cas",        RA(0, 9,10.7), DEC(+1,59, 8,59),  0.525,-0.181, 2.27,"F2III",  REGION_GIANT),
S(98,"Beta Hyi","Hyi",    RA(0,25,45.1), DEC(-1,77,15,15),  2.215, 0.324, 2.80,"G2IV",   REGION_SUBGIANT),
S(99,"Ankaa","Phe",       RA(0,26,17.0), DEC(-1,42,18,22),  0.203,-0.396, 2.39,"K0III",  REGION_GIANT), 
S(168,"Schedar","Cas",     RA(0,40,30.5), DEC(+1,56,32,14),  0.053,-0.032, 2.23,"K0III",  REGION_GIANT),
S(188,"Diphda","Cet",      RA(0,43,35.4), DEC(-1,17,59,12),  0.234, 0.033, 2.04,"G9.5I",  REGION_GIANT),
S(264,"Gamma Cas","Cas",   RA(0,56,42.5), DEC(+1,60,43, 0),  0.026,-0.005, 2.47,"B0IVe",  REGION_SUBGIANT),
S(337,"Mirach","And",      RA(1, 9,43.9), DEC(+1,35,37,14),  0.178,-0.114, 2.06,"M0+II",  REGION_GIANT),
S(424,"Polaris","UMi",     RA(2,31,48.7), DEC(+1,89,15,51),  0.038,-0.015, 2.02,"F7:Ib",  REGION_SUPERGIANT),
S(472,"Achernar","Eri",    RA(1,37,42.9), DEC(-1,57,14,12),  0.095,-0.035, 0.46,"B3Vpe",  REGION_MAIN_SEQUENCE),
S(603,"Almach","And",      RA(2, 3,54.0), DEC(+1,42,19,47),  0.045,-0.052, 2.26,"K3-II",  REGION_GIANT),
S(617,"Hamal","Ari",       RA(2, 7,10.4), DEC(+1,23,27,45),  0.190,-0.148, 2.00,"K2-II",  REGION_GIANT),
S(911,"Menkar","Cet",      RA(3, 2,16.8), DEC(+1, 4, 5,23), -0.009,-0.078, 2.53,"M1.5I",  REGION_GIANT),
S(936,"Algol","Per",       RA(3, 8,10.1), DEC(+1,40,57,20),  0.004,-0.001, 2.12,"B8V",    REGION_MAIN_SEQUENCE),
S(1017,"Mirfak","Per",      RA(3,24,19.4), DEC(+1,49,51,40),  0.024,-0.025, 1.79,"F5Ib",   REGION_SUPERGIANT),
S(1165,"Alcyone","Tau",     RA(3,47,29.1), DEC(+1,24, 6,18),  0.019,-0.046, 2.87,"B7III",  REGION_GIANT),
S(1203,"Zeta Per","Per",    RA(3,54, 7.9), DEC(+1,31,53, 1),  0.006,-0.010, 2.85,"B1Ib",   REGION_SUPERGIANT),
S(1220,"Epsilon Per","Per", RA(3,57,51.2), DEC(+1,40, 0,37),  0.018,-0.026, 2.89,"B0.5V",  REGION_MAIN_SEQUENCE), 
S(1231,"Zaurak","Eri",      RA(3,58, 1.8), DEC(-1,13,30,31),  0.061,-0.111, 2.95,"M0.5I",  REGION_GIANT),

 /* --- Original 3 < V < 4.5: nearby sun-like and high proper motion --- */
S(509,"Tau Ceti","Cet",    RA(1,44, 4.1), DEC(-1,15,56,15), -1.718, 0.856, 3.50,"G8V",    REGION_MAIN_SEQUENCE),
S(1084,"Epsilon Eri","Eri", RA(3,32,55.8), DEC(-1, 9,27,30), -0.974, 0.021, 3.73,"K2V",    REGION_MAIN_SEQUENCE),
S(1136,"Delta Eri","Eri",   RA(3,43,14.9), DEC(-1, 9,45,48), -0.092, 0.745, 3.54,"K0+IV",  REGION_SUBGIANT), 
S(219,"Eta Cas","Cas",     RA(0,49, 6.0), DEC(+1,57,48,57),  1.100,-0.529, 3.44,"F9V",    REGION_MAIN_SEQUENCE),
S(77,"Zeta Tuc","Tuc",    RA(0,20, 4.3), DEC(-1,64,52,29),  1.705, 1.164, 4.23,"F9V",    REGION_MAIN_SEQUENCE),
S(1008,"82 Eri","Eri",      RA(3,19,55.7), DEC(-1,43, 4,11),  3.043, 0.725, 4.27,"G8V",    REGION_MAIN_SEQUENCE),
S(937,"Iota Per","Per",    RA(3, 9, 4.0), DEC(+1,49,36,48),  1.263,-0.091, 4.05,"G0V",    REGION_MAIN_SEQUENCE),
S(321,"Mu Cas","Cas",      RA(1, 8,16.4), DEC(+1,54,55,13),  3.424,-1.596, 5.17,"G5Vb",   REGION_MAIN_SEQUENCE),
S(753,"HR 753","Cet",      RA(2,36, 4.9), DEC(+1, 6,53,13),  1.811, 1.451, 5.82,"K3V",    REGION_MAIN_SEQUENCE), 
S(1010,"Zeta2 Ret","Ret",   RA(3,18,12.9), DEC(-1,62,30,23),  1.336, 0.645, 5.24,"G1V",    REGION_MAIN_SEQUENCE),

 /* --- Original hot stars: O and B, Be stars, Pleiades --- */
S(1228,"Xi Per","Per",      RA(3,58,57.9), DEC(+1,35,47,28),  0.002, 0.000, 4.04,"O7.5I",  REGION_SUPERGIANT),
S(65,"HR 65","Cas",       RA(0,17,43.0), DEC(+1,51,25,59), -0.005,-0.004, 6.14,"O9III",  REGION_GIANT),
S(496,"Phi Per","Per",     RA(1,43,39.6), DEC(+1,50,41,19),  0.026,-0.014, 4.07,"B2Vep",  REGION_MAIN_SEQUENCE),
S(1122,"Delta Per","Per",   RA(3,42,55.5), DEC(+1,47,47,15),  0.028,-0.034, 3.01,"B5III",  REGION_GIANT),
S(1131,"Atik","Per",        RA(3,44,19.1), DEC(+1,32,17,18),  0.007,-0.011, 3.83,"B1III",  REGION_GIANT),
S(1178,"Atlas","Tau",       RA(3,49, 9.7), DEC(+1,24, 3,12),  0.018,-0.047, 3.63,"B8III",  REGION_GIANT),
S(1142,"Electra","Tau",     RA(3,44,52.5), DEC(+1,24, 6,48),  0.019,-0.046, 3.70,"B6III",  REGION_GIANT),
S(1239,"Lambda Tau","Tau",  RA(4, 0,40.8), DEC(+1,12,29,25), -0.006,-0.012, 3.47,"B3V+A",  REGION_MAIN_SEQUENCE),

 /* --- Original A/F stars and subgiants --- */
S(553,"Sheratan","Ari",    RA(1,54,38.4), DEC(+1,20,48,29),  0.096,-0.111, 2.64,"A5V",    REGION_MAIN_SEQUENCE),
S(544,"Mothallah","Tri",   RA(1,53, 4.9), DEC(+1,29,34,44),  0.011,-0.235, 3.41,"F6IV",   REGION_SUBGIANT),
S(897,"Acamar","Eri",      RA(2,58,15.7), DEC(-1,40,18,17), -0.045, 0.019, 3.24,"A4III",  REGION_GIANT),

 /* --- Original supergiants --- */
S(130,"HR 130","Cas",      RA(0,33, 0.0), DEC(+1,62,55,54),  0.003,-0.003, 4.16,"B1Ia",   REGION_SUPERGIANT),
S(1035,"HR 1035","",     RA(3,29, 4.1), DEC(+1,59,56,25), -0.003,-0.004, 4.21,"B9Ia",   REGION_SUPERGIANT),
S(1040,"HR 1040","",     RA(3,29,54.9), DEC(+1,58,52,43),  0.006,-0.002, 4.54,"A0Ia",   REGION_SUPERGIANT),

 /* --- Original cool giants, variables, carbon and S-type (AGB) --- */
S(3,"HR 3","",        RA(0, 5,20.1), DEC(-1, 5,42,27), -0.009, 0.089, 4.61,"K0III",  REGION_GIANT),
S(1208,"Gamma Hyi","Hyi",   RA(3,47,14.3), DEC(-1,74,14,20),  0.047, 0.114, 3.24,"M2III",  REGION_GIANT),
S(681,"Mira","Cet",        RA(2,19,20.7), DEC(-1, 2,58,39), -0.008,-0.237, 3.04,"M7III",  REGION_AGB_CARBON_S),
S(868,"HR 868","Hor",      RA(2,53,52.9), DEC(-1,49,53,25),  0.135, 0.014, 4.00,"M7III",  REGION_AGB_CARBON_S),
S(423,"HR 423","Scl",      RA(1,26,58.1), DEC(-1,32,32,35), -0.011,-0.031, 5.79,"C6II",   REGION_AGB_CARBON_S),
S(977,"HR 977","Hor",      RA(3,12,33.2), DEC(-1,57,19,18),  0.019, 0.009, 5.74,"C5II",   REGION_AGB_CARBON_S),
S(1105,"HR 1105","",     RA(3,42, 9.3), DEC(+1,63,13, 0), -0.020, 0.017, 5.10,"S3.5/",  REGION_AGB_CARBON_S),


S(2491,"Sirius","CMa",     RA(6,45, 8.9), DEC(-1,16,42,58), -0.546,-1.208,-1.46,"A1V",     REGION_MAIN_SEQUENCE),
S(2326,"Canopus","Car",    RA(6,23,57.1), DEC(-1,52,41,44),  0.020, 0.024,-0.72,"F0II",    REGION_GIANT),
S(5459,"Rigil Kentaurus","Cen", RA(14,39,36.5), DEC(-1,60,50, 2),-3.678, 0.481,-0.01,"G2V",  REGION_MAIN_SEQUENCE),
S(5340,"Arcturus","Boo",   RA(14,15,39.7), DEC(+1,19,10,56),-1.093,-1.999,-0.04,"K1.5III", REGION_GIANT),
S(7001,"Vega","Lyr",       RA(18,36,56.3), DEC(+1,38,47, 1), 0.201, 0.282, 0.03,"A0Vvar",  REGION_MAIN_SEQUENCE),
S(1708,"Capella","Aur",    RA(5,16,41.4), DEC(+1,45,59,53),  0.075,-0.427, 0.08,"G5IIIe",  REGION_GIANT),
S(1713,"Rigel","Ori",      RA(5,14,32.3), DEC(-1, 8,12, 6),  0.002,-0.001, 0.12,"B8Ia",    REGION_SUPERGIANT),
S(2943,"Procyon","CMi",    RA(7,39,18.1), DEC(+1, 5,13,30), -0.716,-1.034, 0.38,"F5IV-V",  REGION_SUBGIANT),
S(2061,"Betelgeuse","Ori",  RA(5,55,10.3), DEC(+1, 7,24,25),  0.024, 0.010, 0.50,"M2Iab",   REGION_SUPERGIANT),
S(7557,"Altair","Aql",     RA(19,50,47.0), DEC(+1, 8,52, 6), 0.536, 0.386, 0.77,"A7V",     REGION_MAIN_SEQUENCE),
S(1457,"Aldebaran","Tau",   RA(4,35,55.2), DEC(+1,16,30,33),  0.063,-0.190, 0.85,"K5III",   REGION_GIANT),
S(6134,"Antares","Sco",    RA(16,29,24.4), DEC(-1,26,25,55),-0.010,-0.023, 1.06,"M1.5Iab", REGION_SUPERGIANT),
S(5056,"Spica","Vir",      RA(13,25,11.6), DEC(-1,11, 9,41),-0.043,-0.032, 0.98,"B1III-I", REGION_GIANT),
S(2990,"Pollux","Gem",     RA(7,45,18.9), DEC(+1,28, 1,34), -0.628,-0.046, 1.14,"K0III",   REGION_GIANT),
S(2891,"Castor","Gem",     RA(7,34,36.0), DEC(+1,31,53,18), -0.192,-0.103, 1.58,"A1V",     REGION_MAIN_SEQUENCE),
S(8728,"Fomalhaut","PsA",  RA(22,57,39.0), DEC(-1,29,37,20), 0.329,-0.165, 1.16,"A3V",     REGION_MAIN_SEQUENCE),
S(7924,"Deneb","Cyg",      RA(20,41,25.9), DEC(+1,45,16,49), 0.002, 0.001, 1.25,"A2Iae",   REGION_SUPERGIANT),
S(3982,"Regulus","Leo",    RA(10, 8,22.3), DEC(+1,11,58, 2),-0.248, 0.006, 1.35,"B7V",     REGION_MAIN_SEQUENCE),
S(1790,"Bellatrix","Ori",  RA(5,25, 7.9), DEC(+1, 6,20,59), -0.008,-0.013, 1.64,"B2III",   REGION_GIANT), 
        S(1948,"Alnitak","Ori",    RA(5,40,45.5), DEC(-1, 1,56,34),  0.003, 0.002, 1.74,"O9.7Ib",  REGION_SUPERGIANT),
S(1903,"Alnilam","Ori",    RA(5,36,12.8), DEC(-1, 1,12, 7),  0.001,-0.001, 1.69,"B0Ia",    REGION_SUPERGIANT),
S(1852,"Mintaka","Ori",    RA(5,32, 0.4), DEC(-1, 0,17,57),  0.004,-0.004, 2.25,"B0III",   REGION_GIANT),
S(2004,"Saiph","Ori",      RA(5,47,45.4), DEC(-1, 9,40,11),  0.002,-0.003, 2.06,"B0.5Iav", REGION_SUPERGIANT),
S(4301,"Dubhe","UMa",      RA(11, 3,43.7), DEC(+1,61,45, 3),-0.136,-0.035, 1.81,"K0III",   REGION_GIANT),
S(5054,"Mizar","UMa",      RA(13,23,55.5), DEC(+1,54,55,31), 0.119,-0.026, 2.23,"A2V",     REGION_MAIN_SEQUENCE),
S(5062,"Alcor","UMa",      RA(13,25,13.5), DEC(+1,54,59,17), 0.120,-0.016, 3.99,"A5V",     REGION_MAIN_SEQUENCE),
S(4905,"Alioth","UMa",     RA(12,54, 1.7), DEC(+1,55,57,35), 0.113,-0.010, 1.76,"A0Vp",    REGION_MAIN_SEQUENCE),
S(5191,"Alkaid","UMa",     RA(13,47,32.4), DEC(+1,49,18,48),-0.121,-0.016, 1.85,"B3V",     REGION_MAIN_SEQUENCE),
S(5563,"Kochab","UMi",     RA(14,50,42.3), DEC(+1,74, 9,20),-0.032, 0.011, 2.07,"K4III",   REGION_GIANT),
S(5291,"Thuban","Dra",     RA(14, 4,23.3), DEC(+1,64,22,30),-0.058, 0.000, 3.67,"A0III",   REGION_GIANT),
S(4730,"Acrux","Cru",      RA(12,26,35.9), DEC(-1,63, 5,57),-0.035,-0.012, 0.77,"B0.5IV",  REGION_SUBGIANT),
S(5267,"Hadar","Cen",      RA(14, 3,49.4), DEC(-1,60,22,23),-0.034,-0.025, 0.61,"B1III",   REGION_GIANT),
S(6527,"Shaula","Sco",     RA(17,33,36.5), DEC(-1,37, 6,14),-0.011,-0.027, 1.62,"B2IV",    REGION_SUBGIANT),
S(3748,"Alphard","Hya",    RA(9,27,35.2), DEC(-1, 8,39,31),-0.014, 0.033, 1.99,"K3III",   REGION_GIANT),
S(4534,"Denebola","Leo",   RA(11,49, 3.6), DEC(+1,14,34,19),-0.499,-0.127, 2.14,"A3V",     REGION_MAIN_SEQUENCE),
S(7796,"Sadr","Cyg",       RA(20,22,13.7), DEC(+1,40,15,24), 0.005, 0.000, 2.23,"F8Ib",    REGION_SUPERGIANT),
S(7417,"Albireo","Cyg",    RA(19,30,43.3), DEC(+1,27,57,35),-0.002,-0.001, 3.05,"K3II",    REGION_GIANT),
S(6556,"Rasalhague","Oph", RA(17,34,56.0), DEC(+1,12,33,36), 0.109,-0.222, 2.08,"A5III",   REGION_GIANT),
S(8308,"Enif","Peg",       RA(21,44,11.1), DEC(+1, 9,52,30), 0.023, 0.001, 2.38,"K2Ib",    REGION_SUPERGIANT),
S(8781,"Markab","Peg",     RA(23, 4,45.6), DEC(+1,15,12,19), 0.060,-0.043, 2.49,"B9III",   REGION_GIANT),
S(8775,"Scheat","Peg",     RA(23, 2,57.4), DEC(+1,28, 4,55), 0.191, 0.138, 2.44,"M2III",   REGION_GIANT),
S(4915,"Cor Caroli","CVn", RA(12,56, 1.7), DEC(+1,38,19, 6),-0.231, 0.051, 2.89,"A0p",     REGION_MAIN_SEQUENCE),
S(4210,"Eta Carinae","Car",RA(10,45, 3.6), DEC(-1,59,41, 4),-0.016, 0.012, 6.21,"LBV",     REGION_SUPERGIANT)
};

/* IAU abbreviation -> full name (extended to 42 constellations) */
static const char *constellation_full(const char *abbr)
{
    static const char *tab[][2] = {
        {"And","Andromeda"},{"Ari","Aries"},{"Cas","Cassiopeia"},{"Cet","Cetus"},
        {"Eri","Eridanus"},{"Hor","Horologium"},{"Hyi","Hydrus"},{"Per","Perseus"},
        {"Phe","Phoenix"},{"Ret","Reticulum"},{"Scl","Sculptor"},{"Tau","Taurus"},
        {"Tri","Triangulum"},{"Tuc","Tucana"},{"UMi","Ursa Minor"},{"CMa","Canis Major"},
        {"Car","Carina"},{"Cen","Centaurus"},{"Boo","Boötes"},{"Lyr","Lyra"},
        {"Aur","Auriga"},{"Ori","Orion"},{"CMi","Canis Minor"},{"Aql","Aquila"},
        {"Sco","Scorpius"},{"Vir","Virgo"},{"Gem","Gemini"},{"PsA","Piscis Austrinus"},
        {"Cyg","Cygnus"},{"Leo","Leo"},{"UMa","Ursa Major"},{"Dra","Draco"},
        {"Cru","Crux"},{"Hya","Hydra"},{"Oph","Ophiuchus"},{"Peg","Pegasus"},
        {"CVn","Canes Venatici"}
    };
    for (unsigned i = 0; i < sizeof tab / sizeof tab[0]; i++)
        if (abbr[0] && tab[i][0][0]==abbr[0] && tab[i][0][1]==abbr[1] && tab[i][0][2]==abbr[2])
            return tab[i][1];
    return "unknown";
}

#define SEED_COUNT ((int)(sizeof SEED / sizeof SEED[0]))
