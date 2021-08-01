/* satmars.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/*
 * Original source was taken from
 * ftp://cyrano-se.obspm.fr/pub/3_solar_system/4_satmars/satmars.for
 * 
 * Manual editing was done after conversion:
 *  - removed MAIN_ function
 *  - eliminated calls for formatted input routines and
 *    replaced then with sscanf() calls (to get rid of libf2c)
 * 
 */

#include "f2c.h"

#include <stdio.h>
#include <string.h>

#define SCAN(width, format, variable) do { \
          strncpy(substr, str, width); substr[width] = 0; \
          sscanf(substr, "%" #width #format, variable); \
          str += width; } while (0)

/* Common Block Declarations */

struct ep1_1_ {
    doublereal pm[3], al[3];
};
struct ep1_2_ {
    doublereal dnu, dgam, de, al[3];
};

#define ep1_1 (*(struct ep1_1_ *) &ep1_)
#define ep1_2 (*(struct ep1_2_ *) &ep1_)

struct ep2_1_ {
    doublereal accep;
};

#define ep2_1 (*(struct ep2_1_ *) &ep2_)

struct ep3_1_ {
    doublereal psi0, pipet, alm, dp;
};

#define ep3_1 (*(struct ep3_1_ *) &ep3_)

struct ep4_1_ {
    integer nb[3], nbv[3];
};

#define ep4_1 (*(struct ep4_1_ *) &ep4_)

struct ed1_1_ {
    doublereal dm[3], ald[3];
};
struct ed1_2_ {
    doublereal dnu, dgam, de, al[3];
};

#define ed1_1 (*(struct ed1_1_ *) &ed1_)
#define ed1_2 (*(struct ed1_2_ *) &ed1_)

struct ed2_1_ {
    doublereal acced;
};
struct ed2_2_ {
    doublereal accep;
};

#define ed2_1 (*(struct ed2_1_ *) &ed2_)
#define ed2_2 (*(struct ed2_2_ *) &ed2_)

struct ed4_1_ {
    integer nbd[3], nbvd[3];
};
struct ed4_2_ {
    integer nb[3], nbv[3];
};

#define ed4_1 (*(struct ed4_1_ *) &ed4_)
#define ed4_2 (*(struct ed4_2_ *) &ed4_)

struct {
    doublereal arg[93]	/* was [31][3] */, freq[93]	/* was [31][3] */, cs[
	    93]	/* was [31][3] */, cc[93]	/* was [31][3] */, freq2[93]	
	    /* was [31][3] */, argv[108]	/* was [36][3] */, freqv[108]	
	    /* was [36][3] */, cvs[108]	/* was [36][3] */, cvc[108]	/* 
	    was [36][3] */, frev2[108]	/* was [36][3] */;
} ep5_;

#define ep5_1 ep5_

struct {
    doublereal arg[114]	/* was [38][3] */, freq[114]	/* was [38][3] */, cs[
	    114]	/* was [38][3] */, cc[114]	/* was [38][3] */, 
	    freq2[114]	/* was [38][3] */, argv[81]	/* was [27][3] */, 
	    freqv[81]	/* was [27][3] */, cvs[81]	/* was [27][3] */, 
	    cvc[81]	/* was [27][3] */, frev2[81]	/* was [27][3] */;
} ed5_;

#define ed5_1 ed5_

/* Initialized data */

struct {
    doublereal e_1;
    } ep2_ = { 9.518e-9 };

struct {
    doublereal e_1[4];
    } ep3_ = { 208.5619, 71.0053, 19.373, -1.81103e-7 };

struct {
    integer e_1[6];
    } ep4_ = { 27, 27, 10, 32, 32, 12 };

struct {
    doublereal e_1[6];
    } ed1_ = { -3.33e-5, 1.03e-4, -2.21e-4, 215.2172, 11.2, 224.01 };

struct {
    doublereal e_1;
    } ed2_ = { -3.77e-10 };

struct {
    integer e_1[6];
    } ed4_ = { 37, 38, 25, 27, 27, 16 };

struct {
    doublereal e_1[6];
    } ep1_ = { 4.994e-4, -3.41e-4, 1.46e-4, 171.916, 125.88, 342.91 };


/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__6 = 6;
static integer c__4 = 4;


/*      *********************************************************** */
/*      *     COMPUTATION OF AN EPHEMERIS OF PHOBOS AND DEIMOS    * */
/*      *              FROM ESAPHO AND ESADE THEORIES             * */
/*      *     (MICHELLE CHAPRONT-TOUZE, OBSERVATOIRE DE PARIS)    * */
/*      *********************************************************** */
/*      * Constants from (Chapront-Touze, 1990,                   * */
/*      *         Astron. Astrophys., 240, 159, table 8)          * */
/*      * Phobos' inertial parameters from (Borderies and Yoder,  * */
/*      *         1990, Astron. Astrophys., 233, 235)             * */
/*      *********************************************************** */
/*      *                 List of subroutines                     * */
/*      * EPHOB : computation of position and velocity of Phobos  * */
/*      *         in various reference frames                     * */
/*      * EDEIM : computation of position and velocity of Deimos  * */
/*      *         in various reference frames                     * */
/*      * Subroutines called by EPHOB and EDEIM:                  * */
/*      * CAL   : computation of position and velocity of Phobos  * */
/*      *         in the reference frame of the theories          * */
/*      * CALD  : computation of position and velocity of Deimos  * */
/*      *         in the reference frame of the theories          * */
/*      * REF   : rotation of the reference frame of the theories * */
/*      *         onto the FK5 reference frame                    * */
/*      * REF50 : rotation of the reference frame of the theories * */
/*      *         onto the FK4 reference frame                    * */
/*      * LEC   : loading ESAPHO series, and new constants in     * */
/*      *         the coefficients and in the mean motions        * */
/*      * LECD  : loading ESADE series, and new constants in      * */
/*      *         the coefficients and in the mean motions        * */
/*      *********************************************************** */

/*      This main stands as an example for use of subroutines */
/*      EPHOB and EDEIM. */
/*      It computes, for two Julian dates (DJ2447558.5 and DJ2451515.0), */
/*      the rectangular coordinates of Phobos and Deimos in three */
/*      reference frames (reference frame of the theories, FK4 and FK5) */


/*      **************************************************************** */
/*      * Description of the commons for ESAPHO theory                 * */
/*      * ep1 : dnu (deg/day), dgam, de = corrections to the mean      * */
/*      *       elements nu, gamma, e of ESAPHO theory                 * */
/*      *       lg, h, pi (deg) = mean mean longitude, mean longitudes * */
/*      *       du node and pericentre for Phobos in J2000.0           * */
/*      * ep2 : coefficient of t**2 in the mean longitude of Phobos    * */
/*      * ep3 : parameters fixing Mars' position in degree and         * */
/*      *       correction to Mars' precession of ESAPHO in deg/day    * */
/*      * ep4 : number of terms in the series for position and         * */
/*      *       velocity of Phobos (except for perturbations by        * */
/*      *       Deimos and planets)                                    * */
/*      **************************************************************** */
/*      * Description of the commons for ESADE theory                  * */
/*      * ed1 : dnu (deg/day), dgam, de = corrections to the mean      * */
/*      *       elements nu, gamma, e of ESADE theory                  * */
/*      *       lg, h, pi (deg) = mean mean longitude, mean longitudes * */
/*      *       du node and pericentre for Deimos in J2000.0           * */
/*      * ed2 : coefficient of t**2 in the mean longitude of Deimos    * */
/*      * ed4 : number of terms in the series for position and         * */
/*      *       velocity of Deimos                                     * */
/*      **************************************************************** */









/* Subroutine */ int ephob_(doublereal *dj, doublereal *x, doublereal *v, 
	integer *ir)
{
    /* Initialized data */

    static integer kle = 1;
    static doublereal d2000 = 2451545.;
    static doublereal day = 86400.;

    static doublereal c__[9]	/* was [3][3] */;
    static integer i__, k;
    static doublereal t, w, dc[9]	/* was [3][3] */, vp[3], wp, xp[3];
    extern /* Subroutine */ int cal_(doublereal *, doublereal *, doublereal *)
	    , lec_(void), ref_(doublereal *, doublereal *, doublereal *), 
	    ref50_(doublereal *, doublereal *, doublereal *);

/*      ************************************************************** */
/*      * Computation of position and velocity of Phobos in various  * */
/*      *          reference frames                                  * */
/*      * Input :  dj = Julian ephemeris date                        * */
/*      *          ir = 1 for the reference frame of the theories,   * */
/*      *             = 2 for the FK4 reference frame,               * */
/*      *             = 3 for the FK5 reference frame                * */
/*      * Output : position x(3) and velocity v(3) of Phobos in the  * */
/*      *          reference frame fixed by ir, units = km and  km/s * */
/*      ************************************************************** */



    /* Parameter adjustments */
    --v;
    --x;

    /* Function Body */

    if (kle == 1) {
	lec_();
	kle = 2;
    }

    t = *dj - d2000;
    cal_(&t, xp, vp);

    if (*ir == 1) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    x[i__] = xp[i__ - 1];
	    v[i__] = vp[i__ - 1] / day;
	}
	return 0;
    }

    if (*ir == 2) {
	ref50_(&t, c__, dc);
    } else {
	ref_(&t, c__, dc);
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	w = 0.;
	wp = 0.;
	for (k = 1; k <= 3; ++k) {
	    w += c__[i__ + k * 3 - 4] * xp[k - 1];
	    wp = wp + c__[i__ + k * 3 - 4] * vp[k - 1] + dc[i__ + k * 3 - 4] *
		     xp[k - 1];
	}
	x[i__] = w;
	v[i__] = wp / day;
    }

    return 0;
} /* ephob_ */




/* Subroutine */ int edeim_(doublereal *dj, doublereal *x, doublereal *v, 
	integer *ir)
{
    /* Initialized data */

    static integer kle = 1;
    static doublereal d2000 = 2451545.;
    static doublereal day = 86400.;

    static doublereal c__[9]	/* was [3][3] */;
    static integer i__, k;
    static doublereal t, w, dc[9]	/* was [3][3] */, vp[3], wp, xp[3];
    extern /* Subroutine */ int ref_(doublereal *, doublereal *, doublereal *)
	    , cald_(doublereal *, doublereal *, doublereal *), lecd_(void), 
	    ref50_(doublereal *, doublereal *, doublereal *);

/*      ************************************************************** */
/*      * Computation of position and velocity of Deimos in various  * */
/*      *          reference frames                                  * */
/*      * Input :  dj = Julian ephemeris date                        * */
/*      *          ir = 1 for the reference frame of the theories,   * */
/*      *             = 2 for the FK4 reference frame,               * */
/*      *             = 3 for the FK5 reference frame                * */
/*      * Output : position x(3) and velocity v(3) of Deimos in the  * */
/*      *          reference frame fixed by ir, units = km and km/s  * */
/*      ************************************************************** */



    /* Parameter adjustments */
    --v;
    --x;

    /* Function Body */

    if (kle == 1) {
	lecd_();
	kle = 2;
    }

    t = *dj - d2000;
    cald_(&t, xp, vp);

    if (*ir == 1) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    x[i__] = xp[i__ - 1];
	    v[i__] = vp[i__ - 1] / day;
	}
	return 0;
    }

    if (*ir == 2) {
	ref50_(&t, c__, dc);
    } else {
	ref_(&t, c__, dc);
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	w = 0.;
	wp = 0.;
	for (k = 1; k <= 3; ++k) {
	    w += c__[i__ + k * 3 - 4] * xp[k - 1];
	    wp = wp + c__[i__ + k * 3 - 4] * vp[k - 1] + dc[i__ + k * 3 - 4] *
		     xp[k - 1];
	}
	x[i__] = w;
	v[i__] = wp / day;
    }

    return 0;
} /* edeim_ */




/* Subroutine */ int cal_(doublereal *t, doublereal *x, doublereal *v)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal a;
    static integer i__, j;
    static doublereal t2, ca, sa, xw, deg;
    static integer jmax;

/*      ************************************************************* */
/*      * Computation of position and velocity of Phobos            * */
/*      *          in the reference frame of the theories           * */
/*      * Input  : t = ephemeris time in days reckoned from J2000.0 * */
/*      * Output : position x(3) in km and velocity v(3) in km/day  * */
/*      ************************************************************* */





    /* Parameter adjustments */
    --v;
    --x;

    /* Function Body */
    deg = .017453292519943292;
    t2 = *t * *t * ep2_1.accep * deg;

    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ep4_1.nb[i__ - 1];
	if (i__ != 3) {
	    jmax += 4;
	}
	xw = 0.;
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {
	    a = ep5_1.arg[j + i__ * 31 - 32] + ep5_1.freq[j + i__ * 31 - 32] *
		     *t + ep5_1.freq2[j + i__ * 31 - 32] * t2;
	    sa = sin(a);
	    ca = cos(a);
	    xw = xw + ep5_1.cs[j + i__ * 31 - 32] * sa + ep5_1.cc[j + i__ * 
		    31 - 32] * ca;
	}
	x[i__] = xw;
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ep4_1.nbv[i__ - 1];
	if (i__ != 3) {
	    jmax += 4;
	}
	xw = 0.;
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {
	    a = ep5_1.argv[j + i__ * 36 - 37] + ep5_1.freqv[j + i__ * 36 - 37]
		     * *t + ep5_1.frev2[j + i__ * 36 - 37] * t2;
	    sa = sin(a);
	    ca = cos(a);
	    xw = xw + ep5_1.cvs[j + i__ * 36 - 37] * sa + ep5_1.cvc[j + i__ * 
		    36 - 37] * ca;
	}
	v[i__] = xw;
    }

    return 0;
} /* cal_ */




/* Subroutine */ int cald_(doublereal *t, doublereal *x, doublereal *v)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal a;
    static integer i__, j;
    static doublereal t2, ca, sa, xw, deg;
    static integer jmax;

/*      ************************************************************* */
/*      * Computation of position and velocity of Deimos            * */
/*      *          in the reference frame of the theories           * */
/*      * Input :  t = ephemeris time in days reckoned from J2000.0 * */
/*      * Output : position x(3) in km and velocity v(3) in km/day  * */
/*      ************************************************************* */





    /* Parameter adjustments */
    --v;
    --x;

    /* Function Body */
    deg = .017453292519943292;
    t2 = *t * *t * ed2_2.accep * deg;

    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ed4_2.nb[i__ - 1];
	xw = 0.;
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {
	    a = ed5_1.arg[j + i__ * 38 - 39] + ed5_1.freq[j + i__ * 38 - 39] *
		     *t + ed5_1.freq2[j + i__ * 38 - 39] * t2;
	    sa = sin(a);
	    ca = cos(a);
	    xw = xw + ed5_1.cs[j + i__ * 38 - 39] * sa + ed5_1.cc[j + i__ * 
		    38 - 39] * ca;
	}
	x[i__] = xw;
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ed4_2.nbv[i__ - 1];
	xw = 0.;
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {
	    a = ed5_1.argv[j + i__ * 27 - 28] + ed5_1.freqv[j + i__ * 27 - 28]
		     * *t + ed5_1.frev2[j + i__ * 27 - 28] * t2;
	    sa = sin(a);
	    ca = cos(a);
	    xw = xw + ed5_1.cvs[j + i__ * 27 - 28] * sa + ed5_1.cvc[j + i__ * 
		    27 - 28] * ca;
	}
	v[i__] = xw;
    }

    return 0;
} /* cald_ */




/* Subroutine */ int ref_(doublereal *t, doublereal *c__, doublereal *dc)
{
    /* Initialized data */

    static doublereal g[9]	/* was [3][3] */ = { 1.,-4.79966e-7,0.,
	    4.4036e-7,.917482137087,.397776982902,-1.90919e-7,-.397776982902,
	    .917482137087 };

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal d__[9]	/* was [3][3] */;
    static integer i__, j, k;
    static doublereal q, w, ah, am[9]	/* was [3][3] */, ar[9]	/* was [3][3] 
	    */, dp[9]	/* was [3][3] */, tm, qp, wp, gm2, deg, gam, rac, rad,
	     ahp, amp[9]	/* was [3][3] */, ggp, csh, arp[9]	/* 
	    was [3][3] */, csq, snh, cst, thp, tet, snq, snt, coef, gamp, 
	    csth, tetp, snth, coefp;

/*      ************************************************************* */
/*      * Rotation of the reference frame of the theories onto      * */
/*      *          the FK5 reference frame                          * */
/*      * Input  : t = ephemeris time in days reckoned from J2000.0 * */
/*      * Output : c(3,3)  = rotation matrix                        * */
/*      *          dc(3,3) = derivative of the rotation matrix      * */
/*      ************************************************************* */



    /* Parameter adjustments */
    dc -= 4;
    c__ -= 4;

    /* Function Body */

    tm = *t / 365250.;
    rad = 4.848136811095359e-6;
    deg = rad * 3600;

    ah = rad * 178409.13618 + tm * (tm * (-.001117775392103901 - tm * 
	    3.427183852553758e-5) - .05149158068948755);
    gam = tm * (tm * (-1.968131454406586e-5 - tm * 2.505007528551377e-7) - 
	    7.10928404247926e-4) + .01614120767052974;
    tet = (*t * 2.507593e-6 + 35.496817571) * deg;

    ahp = (tm * (-.0022355507842078022 - tm * 1.0281551557661273e-4) - 
	    .05149158068948755) / 365250.;
    gamp = (tm * (-3.9362629088131719e-5 - tm * 7.5150225856541298e-7) - 
	    7.10928404247926e-4) / 365250.;
    tetp = deg * 2.507593e-6;
    thp = ahp + tetp;
    ggp = gam * 4 * gamp;

    csh = cos(ah);
    snh = sin(ah);
    cst = cos(tet);
    snt = sin(tet);
    csth = csh * cst - snh * snt;
    snth = snt * csh + snh * cst;
    gm2 = gam * gam;
    rac = sqrt(1 - gm2);
    gm2 *= 2;

    am[0] = csth + gm2 * snh * snt;
    am[3] = -snth + gm2 * snh * cst;
    am[1] = snth - gm2 * csh * snt;
    am[4] = csth - gm2 * csh * cst;
    coef = gam * 2 * rac;
    am[6] = coef * snh;
    am[7] = -coef * csh;
    am[2] = coef * snt;
    am[5] = coef * cst;
    am[8] = 1 - gm2;

    amp[0] = -thp * snth + ggp * snh * snt + gm2 * (ahp * csh * snt + tetp * 
	    snh * cst);
    amp[3] = -thp * csth + ggp * snh * cst + gm2 * (ahp * csh * cst - tetp * 
	    snh * snt);
    amp[1] = thp * csth - ggp * csh * snt + gm2 * (ahp * snh * snt - tetp * 
	    csh * cst);
    amp[4] = -thp * snth - ggp * csh * cst + gm2 * (ahp * snh * cst + tetp * 
	    csh * snt);
    coefp = gamp * 2 * (1 - gm2) / rac;
    amp[6] = coef * ahp * csh + coefp * snh;
    amp[7] = coef * ahp * snh - coefp * csh;
    amp[2] = coef * tetp * cst + coefp * snt;
    amp[5] = -coef * tetp * snt + coefp * cst;
    amp[8] = gam * -4 * gamp;

    q = (*t * 3.269878e-7 + 25.19202802) * deg;
    qp = deg * 3.269878e-7;
    csq = cos(q);
    snq = sin(q);

    ar[0] = 1.;
    for (i__ = 2; i__ <= 3; ++i__) {
	ar[i__ * 3 - 3] = 0.;
	ar[i__ - 1] = 0.;
    }
    ar[4] = csq;
    ar[7] = -snq;
    ar[5] = snq;
    ar[8] = csq;

    for (i__ = 1; i__ <= 3; ++i__) {
	arp[i__ * 3 - 3] = 0.;
	arp[i__ - 1] = 0.;
    }
    arp[4] = -qp * snq;
    arp[7] = -qp * csq;
    arp[5] = qp * csq;
    arp[8] = -qp * snq;

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    w = 0.;
	    wp = 0.;
	    for (k = 1; k <= 3; ++k) {
		wp = wp + amp[i__ + k * 3 - 4] * ar[k + j * 3 - 4] + am[i__ + 
			k * 3 - 4] * arp[k + j * 3 - 4];
		w += am[i__ + k * 3 - 4] * ar[k + j * 3 - 4];
	    }
	    dp[i__ + j * 3 - 4] = wp;
	    d__[i__ + j * 3 - 4] = w;
	}
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    w = 0.;
	    wp = 0.;
	    for (k = 1; k <= 3; ++k) {
		wp += g[i__ + k * 3 - 4] * dp[k + j * 3 - 4];
		w += g[i__ + k * 3 - 4] * d__[k + j * 3 - 4];
	    }
	    dc[i__ + j * 3] = wp;
	    c__[i__ + j * 3] = w;
	}
    }

    return 0;
} /* ref_ */




/* Subroutine */ int ref50_(doublereal *t, doublereal *c__, doublereal *dc)
{
    /* Initialized data */

    static doublereal g[9]	/* was [3][3] */ = { .999925674124,
	    -.011181963465,-.004859004081,.01219205172,.917413967951,
	    .39774736364,1.0121726e-5,-.397777041948,.917482111431 };

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal d__[9]	/* was [3][3] */;
    static integer i__, j, k;
    static doublereal q, w, ah, am[9]	/* was [3][3] */, ar[9]	/* was [3][3] 
	    */, dp[9]	/* was [3][3] */, tm, qp, wp, gm2, deg, gam, rac, rad,
	     ahp, amp[9]	/* was [3][3] */, ggp, csh, arp[9]	/* 
	    was [3][3] */, csq, snh, cst, thp, tet, snq, snt, coef, gamp, 
	    csth, tetp, snth, coefp;

/*      ************************************************************* */
/*      * Rotation of the reference frame of the theories onto      * */
/*      *          the FK4 reference frame                          * */
/*      * Input  : t = ephemeris time in days reckoned from J2000.0 * */
/*      * Output : c(3,3)  = rotation matrix                        * */
/*      *          dc(3,3) = derivative of the rotation matrix      * */
/*      ************************************************************* */



    /* Parameter adjustments */
    dc -= 4;
    c__ -= 4;

    /* Function Body */
    tm = *t / 365250.;
    rad = 4.848136811095359e-6;
    deg = rad * 3600;

    ah = rad * 178409.13618 + tm * (tm * (-.001117775392103901 - tm * 
	    3.427183852553758e-5) - .05149158068948755);
    gam = tm * (tm * (-1.968131454406586e-5 - tm * 2.505007528551377e-7) - 
	    7.10928404247926e-4) + .01614120767052974;
    tet = (*t * 2.507593e-6 + 35.496817571) * deg;

    ahp = (tm * (-.0022355507842078022 - tm * 1.0281551557661273e-4) - 
	    .05149158068948755) / 365250.;
    gamp = (tm * (-3.9362629088131719e-5 - tm * 7.5150225856541298e-7) - 
	    7.10928404247926e-4) / 365250.;
    tetp = deg * 2.507593e-6;
    thp = ahp + tetp;
    ggp = gam * 4 * gamp;

    csh = cos(ah);
    snh = sin(ah);
    cst = cos(tet);
    snt = sin(tet);
    csth = csh * cst - snh * snt;
    snth = snt * csh + snh * cst;
    gm2 = gam * gam;
    rac = sqrt(1 - gm2);
    gm2 *= 2;

    am[0] = csth + gm2 * snh * snt;
    am[3] = -snth + gm2 * snh * cst;
    am[1] = snth - gm2 * csh * snt;
    am[4] = csth - gm2 * csh * cst;
    coef = gam * 2 * rac;
    am[6] = coef * snh;
    am[7] = -coef * csh;
    am[2] = coef * snt;
    am[5] = coef * cst;
    am[8] = 1 - gm2;

    amp[0] = -thp * snth + ggp * snh * snt + gm2 * (ahp * csh * snt + tetp * 
	    snh * cst);
    amp[3] = -thp * csth + ggp * snh * cst + gm2 * (ahp * csh * cst - tetp * 
	    snh * snt);
    amp[1] = thp * csth - ggp * csh * snt + gm2 * (ahp * snh * snt - tetp * 
	    csh * cst);
    amp[4] = -thp * snth - ggp * csh * cst + gm2 * (ahp * snh * cst + tetp * 
	    csh * snt);
    coefp = gamp * 2 * (1 - gm2) / rac;
    amp[6] = coef * ahp * csh + coefp * snh;
    amp[7] = coef * ahp * snh - coefp * csh;
    amp[2] = coef * tetp * cst + coefp * snt;
    amp[5] = -coef * tetp * snt + coefp * cst;
    amp[8] = gam * -4 * gamp;

    q = (*t * 3.269878e-7 + 25.19202802) * deg;
    qp = deg * 3.269878e-7;
    csq = cos(q);
    snq = sin(q);

    ar[0] = 1.;
    for (i__ = 2; i__ <= 3; ++i__) {
	ar[i__ * 3 - 3] = 0.;
	ar[i__ - 1] = 0.;
    }
    ar[4] = csq;
    ar[7] = -snq;
    ar[5] = snq;
    ar[8] = csq;

    for (i__ = 1; i__ <= 3; ++i__) {
	arp[i__ * 3 - 3] = 0.;
	arp[i__ - 1] = 0.;
    }
    arp[4] = -qp * snq;
    arp[7] = -qp * csq;
    arp[5] = qp * csq;
    arp[8] = -qp * snq;

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    w = 0.;
	    wp = 0.;
	    for (k = 1; k <= 3; ++k) {
		wp = wp + amp[i__ + k * 3 - 4] * ar[k + j * 3 - 4] + am[i__ + 
			k * 3 - 4] * arp[k + j * 3 - 4];
		w += am[i__ + k * 3 - 4] * ar[k + j * 3 - 4];
	    }
	    dp[i__ + j * 3 - 4] = wp;
	    d__[i__ + j * 3 - 4] = w;
	}
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    w = 0.;
	    wp = 0.;
	    for (k = 1; k <= 3; ++k) {
		wp += g[i__ + k * 3 - 4] * dp[k + j * 3 - 4];
		w += g[i__ + k * 3 - 4] * d__[k + j * 3 - 4];
	    }
	    dc[i__ + j * 3] = wp;
	    c__[i__ + j * 3] = w;
	}
    }

    return 0;
} /* ref50_ */




/* Subroutine */ int lec_(void)
{
    /* Initialized data */

    static doublereal asup[4] = { 183.113,199.2812,255.194,271.3608 };
    static doublereal fsup[4] = { 1128.826758,-1128.862761,1128.849665,
	    -1128.839842 };
    static doublereal rx1s[4] = { 51.71,-51.71,63.305,-63.208 };
    static doublereal rx2s[4] = { 51.71,51.71,63.305,63.209 };
    static doublereal fipsi = 350.8919885;
    static doublereal fipi = 1.772311e-5;
    static doublereal fid = 1128.320721;
    static doublereal fif = 1129.280784;
    static doublereal fil = 1128.409439;
    static doublereal filp = .5240207;
    static doublereal fnu = 1128.84426;
    static struct {
	char e_1[4636];
	char fill_2[1220];
	} equiv_18 = { "  0  1  1  0  0  1-184667147.923 -61587324  3570444 "
		"  2774147  0  1  1  0 -1  1      1606.495      2670     -278"
		"    107003  0  1  1  0  1  1  -2771993.617   -927369    5389"
		"0-184737001  0  1  1 -2  0  1     17230.055      5772  35680"
		"97      -255  0  1  1  0  2  1    -46798.873    -15701      "
		"884  -6238465  0  1  1  0  0  0     -8832.492      6011     "
		" 177       157  0  1  1  0  0  2      8838.522     -6027    "
		" -197      -118  0  2  2 -1  0  3     22776.763    163962   "
		"649508      1773  0  1  3  0 -1  1    -12763.760   -136428  "
		"   9848   -842940  0  0  0  1  0 -1    -22555.148   -161818 "
		" -626771     -1714  0  1 -1  0  0  1      5700.884     -3171"
		"      -19      6279  0  3  1  0  0  3      5697.288     -377"
		"3     -591       -78  0  1  1  0 -2  1      5221.378      17"
		"70      -97    696244  0  1  1  1  0  1         0.229       "
		"  0        0         0  0  2  2 -1  0  2      2846.423     -"
		"8794   301896      -197  0  0  0  1  0  0     -2846.188     "
		" 8820  -301833       157  2 -3 -3  0  0 -3      5485.491    "
		"  2702     -281      -140  2 -1 -1  0  0 -1      -667.307   "
		"   1195      193         7  1  0  0  0  0  0      -532.706  "
		"   -1880      153       -79  0  1  3  0 -2  1      2677.804 "
		"    39842    -1221    212926  0  1 -1  0  2  1     -2486.738"
		"    -37783     1083   -187395  0  3  1  0  0  4      1248.00"
		"8      -785        0         0  0  1 -1  0  0  2      1237.2"
		"87      -748      118       511  1 -2 -2  0  0 -2     -1721."
		"139     -8128      499        66  0  1  1  0  0 -1      -617"
		".237       405        0         0  0  1  1  0  0  3       61"
		"6.876      -428        0         0  0  3  1  0  0  2      -5"
		"43.548       343        0        39  0  1 -1  0  0  0      -"
		"540.379       345      -19       177  0  1  3  0 -1  0     -"
		"1054.707     -2316      196    -70261  0  2  0  1  0  2     "
		"  517.078      -410    53482         0  3 -4 -4  0  0 -4    "
		"  -845.723     -5333       60       302  0  1  1  0  3  1   "
		"   -832.288      -280        0   -166394  0  1  1  0  0  1 1"
		"84667155.134  61587311 -3570720  -2773871  0  1  1  0 -1  1 "
		"    -1606.803     -2670      278   -107023  0  1  1  0  1  1"
		"   2771993.601    927369   -53890 184737001  0  1  1 -2  0  "
		"1    -17237.270     -5759 -3569615       255  0  1  1  0  2 "
		" 1     46798.873     15701     -884   6238465  0  1  1  0  0"
		"  0      8834.010     -6030     -177       -78  0  1  1  0  "
		"0  2     -8834.307      6051       98       118  0  2  2 -1 "
		" 0  3    -22776.766   -163962  -649508     -1773  0  1  3  0"
		" -1  1     12763.760    136428    -9848    842940  0  0  0  "
		"1  0 -1     22592.199    161756   630613      1714  0  3  1 "
		" 0  0  3     -5696.686      3774      611        78  0  1 -1"
		"  0  0  1      5683.941     -4349    -1161     -6456  0  1  "
		"1  1  0  1        -0.229         0        0         0  0  1 "
		" 1  0 -2  1     -5197.382     -1766       97   -693055  0  0"
		"  0  1  0  0      3410.371    -10581   362440      -275  0  "
		"2  2 -1  0  2     -2846.415      8794  -301896       197  2 "
		"-3 -3  0  0 -3      5485.491      2702     -281      -140  2"
		" -1 -1  0  0 -1       667.306     -1195     -193        -7  "
		"1  0  0  0  0  0       532.706      1880     -153        79 "
		" 0  1  3  0 -2  1     -2677.794    -39842     1221   -212926"
		"  0  1 -1  0  2  1      2487.540     37784    -1083    18751"
		"3  0  1 -1  0  0  2      1255.765      -827     -157      -5"
		"11  0  3  1  0  0  4     -1248.014       785        0       "
		"  0  1 -2 -2  0  0 -2     -1721.139     -8128      499      "
		"  66  0  1  1  0  0 -1       617.168      -405        0     "
		"    0  0  1  1  0  0  3      -617.663       408        0    "
		"     0  0  1 -1  0  0  0      -547.417       358       39   "
		"   -118  0  3  1  0  0  2       543.495      -343        0  "
		"     -39  0  2  0  1  0  2      -538.076       389   -55652 "
		"        0  0  1  3  0 -1  0      1054.707      2316     -196"
		"     70261  3 -4 -4  0  0 -4      -845.723     -5333       6"
		"0       302  0  1  1  0  3  1       832.288       280       "
		" 0    166394  0  0  0  1  0  0   3567800.464   1190831369388"
		"343    -53586  0  1  1  0  0  1    -29184.285     91687   -1"
		"2826      2246  0  0  0  1  1  0     53599.765     18017  55"
		"49359   3572140  0  1  1  0  0  2     20782.155    221432   "
		"-29113      5775  0  1 -1  0  0  1      8999.931    -12132  "
		"   1850      -275  0  0  0  0  1  0         0.146         0 "
		"       0         0  0  1 -1  0  0  2      2181.635     -2517"
		"      373       -39  0  1  1  0  0  0     -1913.567      330"
		"4     -472        78  0  3  1  0  0  3      1089.188      10"
		"88      -78         0  2 -2 -2  0  0 -2        13.355       "
		" 19        0         0  0  1 -1  0  0  0      -664.567      "
		"1146     -177        19  0  0  0  1  2  0       905.431     "
		"  314    93730    120442" };

    static struct {
	char e_1[4256];
	char fill_2[1120];
	} equiv_19 = { "         -0.002     0       0      0     0          "
		"1103          0.000     0       0      0     0           -37"
		"          0.000     0       0      0     0           -12    "
		"      0.000     0       0      0     0             0        "
		"  0.000     0       0      0     0             0          0."
		"022     0       0      0     0           472          0.022 "
		"    0       0      0     0          -453          2.366     "
		"0       0      0     0          1536         -0.035     0   "
		"    0      0     0          -118          2.362     0       "
		"0      0     0         -1516         -0.002     0       0   "
		"   0     0          -295         -0.002     0       0      0"
		"     0          -295          0.000     0       0      0    "
		" 0             0     -10422.818-10386-1007363    630     0  "
		"           0         -0.156     0       0      0     0      "
		"  136994         -0.205     0       0      0     0       -13"
		"7041          2.226     2       0      0     0             0"
		"        -47.485  -125       0      0     0             0    "
		"    -52.948  -134       6     -6     0             0        "
		"  0.203     0       0      0     0            19          0."
		"203     0       0      0     0           -19          0.000 "
		"    0       0      0     0           -39          0.000     "
		"0       0      0     0           -39       -165.614  -628   "
		"   33      0     0             0          0.000     0       "
		"0      0     0             0          0.000     0       0   "
		"   0     0             0          0.000     0       0      0"
		"     0            59          0.000     0       0      0    "
		" 0            59          0.000     0       0      0     0  "
		"           0          0.000     0       0      0     0      "
		"       0        135.835   841       0      0     0          "
		"   0          0.000     0       0      0     0             0"
		"         -0.002     0       0      0     0         -1438    "
		"      0.000     0       0      0     0            37        "
		"  0.000     0       0      0     0            12          0."
		"000     0       0      0     0             0          0.000 "
		"    0       0      0     0             0          0.022     "
		"0       0      0     0          -492          0.022     0   "
		"    0      0     0           689          2.366     0       "
		"0      0     0         -1536         -0.035     0       0   "
		"   0     0           118          2.362     0       0      0"
		"     0          1516         -0.002     0       0      0    "
		" 0           315          0.000     0       0      0     0  "
		"        -196     -10422.818-10386-1007363    630     0      "
		"       0          0.000     0       0      0     0          "
		"   0         -0.156     0       0      0     0        176342"
		"         -0.156     0       0      0     0       -136994    "
		"     -2.226    -2       0      0     0             0        "
		"-47.485  -125       0      0     0             0        -52."
		"948  -134       6     -6     0             0          0.203 "
		"    0       0      0     0           -19          0.203     "
		"0       0      0     0            19          0.000     0   "
		"    0      0     0           -19          0.000     0       "
		"0      0     0            39        165.614   628     -33   "
		"   0     0             0          0.000     0       0      0"
		"     0             0          0.000     0       0      0    "
		" 0             0          0.000     0       0      0     0  "
		"          39          0.000     0       0      0     0      "
		"     -59          0.000     0       0      0     0          "
		"   0          0.000     0       0      0     0             0"
		"       -135.835  -841       0      0     0             0    "
		"      0.000     0       0      0     0             0        "
		"  0.000     0       0      0     0           -26         -2."
		"549    -2       0      0     0      -2035142          0.000 "
		"    0       0      0     0             0         -0.363     "
		"0       0      0     0          2089          0.000     0   "
		"    0      0     0           -98      -8100.500 -8096     78"
		"7-504217     0             0          0.000     0       0   "
		"   0     0           -19          0.000     0       0      0"
		"     0           -39          0.000     0       0      0    "
		" 0            39      -1306.279 -1135     271   -108     0  "
		"           0          0.000     0       0      0     0      "
		"       0          0.000     0       0      0     0          "
		"   0" };

    static struct {
	char e_1[3648];
	char fill_2[969];
	} equiv_8 = { "  0  1  1  0  0  1 9372991.9756-6247057  -181222  -14"
		"0805  0  1  1  0 -1  1 -210991.0229  140430     4098-1406605"
		"9  0  1  1  0  1  1   70361.4467  -46804    -1370  4689174  "
		"0  1 -1  0  1  1    -970.2930   -9509      752   -64574  0  "
		"1  1 -2  0  1     873.8556    -582   180963      -13  0  1  "
		"1 -1  0  1       0.0000       0        0        0  0  1  1  "
		"0  2  1     791.9808    -526      -15   105574  0  1  1  0  "
		"0  0     448.5114    -754       -9       -8  0  1  1  0  0  "
		"2    -448.4011     754       10        6  0  2  2 -1  0  3  "
		" -1155.9707   -7995   -32964      -90  0  1  3  0 -1  1     "
		"324.0074    3139     -250    21398  0  0  0  1  0 -1    1144"
		".9013    7898    31815       87  0  1 -1  0  0  1     289.62"
		"38    -451       -1      319  0  3  1  0  0  3    -288.9041 "
		"    480       30        4  0  1  1  0 -2  1     265.2220    "
		"-175       -5    35366  0  1  1  1  0  1      -0.0058       "
		"0        0        0  0  2  2 -1  0  2    -144.5293     591  "
		" -15329       10  0  0  0  1  0  0     144.4057    -592    1"
		"5314       -8  2 -3 -3  0  0 -3     117.0670     -90       -"
		"6       -3  2 -1 -1  0  0 -1     -89.5280     397       26  "
		"      1  1  0  0  0  0  0      86.9834     307      -25     "
		"  13  0  1  3  0 -2  1    -135.9361   -1940       62   -1080"
		"9  0  1 -1  0  1  2     -79.1722     -99       16    -5280  "
		"0  1 -1  0  2  1     126.1976    1845      -55     9510  0  "
		"3  1  0  0  4     -63.2560     103        0        0  0  1 -"
		"1  0  0  2      62.8875    -101        6       26  1 -2 -2  "
		"0  0 -2     -51.7171    -183       15        2  0  1  1  0  "
		"0  1 9372992.3416-6247058  -181236  -140791  0  1  1  0 -1  "
		"1 -211031.5150  140592     4090-14068768  0  1  1  0  1  1  "
		" 70361.4463  -46804    -1370  4689174  0  1 -1  0  1  1    -"
		"980.1755   -9494      752   -65240  0  1  1 -2  0  1     874"
		".2215    -583   181040      -13  0  1  1 -1  0  1       0.00"
		"00       0        0        0  0  1  1  0  2  1     791.9808 "
		"   -526      -15   105574  0  1  1  0  0  0     448.5885    "
		"-755       -9       -4  0  1  1  0  0  2    -448.1873     75"
		"5        5        6  0  2  2 -1  0  3   -1155.9709   -7995  "
		" -32964      -90  0  1  3  0 -1  1     324.0074    3139     "
		"-250    21398  0  0  0  1  0 -1    1146.7820    7893    3201"
		"0       87  0  3  1  0  0  3    -288.8736     480       31  "
		"      4  0  1 -1  0  0  1    -288.7630     510       59     "
		" 328  0  1  1  1  0  1      -0.0058       0        0        "
		"0  0  1  1  0 -2  1     264.0031    -174       -5    35204  "
		"0  0  0  1  0  0     173.0304    -710    18389      -14  0  "
		"2  2 -1  0  2    -144.5289     591   -15329       10  2 -3 -"
		"3  0  0 -3    -117.0670      90        6        3  2 -1 -1  "
		"0  0 -1     -89.5279     397       26        1  1  0  0  0  "
		"0  0      86.9835     307      -25       13  0  1  3  0 -2  "
		"1    -135.9356   -1940       62   -10809  0  1 -1  0  1  2  "
		"   -81.3238     -96       16    -5426  0  1 -1  0  2  1     "
		"126.2383    1845      -55     9516  0  1 -1  0  0  2     -63"
		".8267     106        8       26  0  3  1  0  0  4     -63.25"
		"63     103        0        0  1 -2 -2  0  0 -2      51.7171 "
		"    183      -15       -2  0  0  0  1  0  0  181017.8207 -12"
		"0692 18741486    -2723  0  0  0  1 -1  0   -4077.2629    271"
		"3  -422134  -271818  0  1  1  0  0  1   -1481.2817    6135  "
		"   -651      114  0  0  0  1  1  0    1360.2586    -903   14"
		"0832    90654  0  1  1  0  0  2    1054.3326   10180    -147"
		"7      293  0  1 -1  0  0  1    -457.2263    1074      -94  "
		"     14  0  0  0  0  1  0       0.0074       0        0     "
		"   0  0  1 -1  0  0  2    -110.8858     239      -19        "
		"2  0  1  1  0  0  0     -97.1704     265      -24        4  "
		"0  3  1  0  0  3      55.2317       0       -4        0" };

    static struct {
	char e_1[3264];
	char fill_2[867];
	} equiv_9 = { "       -0.0001     0     0     0     0          -56  "
		"     -0.0002     0     0     0     0            1        0.0"
		"000     0     0     0     0            0       -0.0028     0"
		"     0     0     0          -11        0.0000     0     0   "
		"  0     0            0     -794.4492     3-76760   -34     0"
		"            0        0.0000     0     0     0     0         "
		"   0        0.0011     0     0     0     0          -24     "
		"   0.0011     0     0     0     0           23        0.1201"
		"     0     0     0     0          -78       -0.0009     0   "
		"  0     0     0            3        0.1199     0     0     0"
		"     0           77        0.0001     0     0     0     0   "
		"       -15       -0.0001     0     0     0     0           1"
		"5        0.0000     0     0     0     0            0     -26"
		"4.4599     1-25560    16     0            0       -0.0079   "
		"  0     0     0     0        -6956       -0.0104     0     0"
		"     0     0         6953       -0.0475     0     0     0   "
		"  0            0        6.3708     0     0     0     0      "
		"      0       -8.6457   -22     1    -1     0            0  "
		"      0.0103     0     0     0     0           -1        0.0"
		"000     0     0     0     0           -1        0.0103     0"
		"     0     0     0            1        0.0000     0     0   "
		"  0     0            2        0.0000     0     0     0     0"
		"           -2        4.9764    13    -1     0     0         "
		"   0        0.0001     0     0     0     0          -73     "
		"  -0.0001     0     0     0     0            2        0.0000"
		"     0     0     0     0            0        0.0029     0   "
		"  0     0     0          -11        0.0000     0     0     0"
		"     0            0      794.4491    -3 76760    34     0   "
		"         0        0.0000     0     0     0     0            "
		"0       -0.0011     0     0     0     0          -25       -"
		"0.0011     0     0     0     0           35       -0.1201   "
		"  0     0     0     0          -78        0.0009     0     0"
		"     0     0            3       -0.1199     0     0     0   "
		"  0           77        0.0001     0     0     0     0      "
		"     16        0.0000     0     0     0     0           10  "
		"    264.4599    -1 25560   -16     0            0        0.0"
		"000     0     0     0     0            0        0.0079     0"
		"     0     0     0         8947        0.0079     0     0   "
		"  0     0        -6956       -0.0475     0     0     0     0"
		"            0       -6.3708     0     0     0     0         "
		"   0        8.6457    22    -1     1     0            0     "
		"  -0.0103     0     0     0     0           -1        0.0000"
		"     0     0     0     0           -1       -0.0103     0   "
		"  0     0     0            1        0.0000     0     0     0"
		"     0            1        0.0000     0     0     0     0   "
		"         2        4.9764    13    -1     0     0            "
		"0        0.0000     0     0     0     0           -3        "
		"0.0000     0     0     0     0            0        0.1294   "
		"  0     0     0     0      -103296        0.0000     0     0"
		"     0     0            0        0.0184     0     0     0   "
		"  0          106        0.0000     0     0     0     0      "
		"      5      411.3091     0   -40 25602     0            0  "
		"      0.0000     0     0     0     0            1        0.0"
		"000     0     0     0     0           -2        0.0000     0"
		"     0     0     0            2" };


    /* System generated locals */
    integer i__1;
    icilist ici__1;


    /* Local variables */
    static doublereal f[6];
    static integer i__, j, k;
#define v ((char *)&equiv_18)
#define w ((char *)&equiv_19)
#define x ((char *)&equiv_8)
#define y ((char *)&equiv_9)
    static doublereal a0[6], r1, r2;
#define v1 ((char *)&equiv_18)
#define v2 ((char *)&equiv_18 + 1952)
#define x1 ((char *)&equiv_8)
#define x2 ((char *)&equiv_8 + 1539)
#define x3 ((char *)&equiv_8 + 3078)
#define y1 ((char *)&equiv_9)
#define y2 ((char *)&equiv_9 + 1377)
#define y3 ((char *)&equiv_9 + 2754)
#define v3 ((char *)&equiv_18 + 3904)
#define w1 ((char *)&equiv_19)
#define w2 ((char *)&equiv_19 + 1792)
#define w3 ((char *)&equiv_19 + 3584)
    static doublereal ai[6], pf[6], aw;
    static integer jp;
    static doublereal fw;
    static integer ic1[4], ic2[4];
#define v1p ((char *)&equiv_18 + 1098)
#define v2p ((char *)&equiv_18 + 3050)
#define x1p ((char *)&equiv_8 + 1026)
#define x2p ((char *)&equiv_8 + 2565)
#define y1p ((char *)&equiv_9 + 918)
#define y2p ((char *)&equiv_9 + 2295)
#define w1p ((char *)&equiv_19 + 1008)
#define w2p ((char *)&equiv_19 + 2800)
    static doublereal deg;
    static integer iar[6];
    static doublereal dpn;
    static integer jmax;
    static doublereal dnun;

/*      ********************************************************* */
/*      * Transformation of the series from ESAPHO.             * */
/*      * Introduction of the new constants in the coefficients * */
/*      * and in the mean motions                               * */
/*      ********************************************************* */







/*      **************************** */
/*      * Coefficients from ESAPHO * */
/*      **************************** */

/*      *************************************** */
/*      * Perturbations by Deimos and planets * */
/*      *************************************** */

/*      ************************ */
/*      * Arguments in J2000.0 * */
/*      ************************ */
    ai[0] = ep3_1.psi0;
    ai[1] = ep3_1.pipet;
    ai[2] = ep1_2.al[0] - ep3_1.pipet - ep3_1.alm;
    ai[3] = ep1_2.al[0] - ep1_2.al[1];
    ai[4] = ep1_2.al[0] - ep1_2.al[2];
    ai[5] = ep3_1.alm;

/*      **************************** */
/*      * Mean motions from ESAPHO * */
/*      **************************** */

    dpn = ep3_1.dp * 1e5 / fnu;
    dnun = ep1_2.dnu / fnu;
    pf[0] = fipsi;
    pf[1] = fipi;
    pf[2] = fid;
    pf[3] = fif;
    pf[4] = fil;
    pf[5] = filp;

/*      ************************* */
/*      * Conversion to radians * */
/*      ************************* */
    deg = .017453292519943292;
    for (k = 1; k <= 6; ++k) {
	f[k - 1] = pf[k - 1] * deg;
	a0[k - 1] = ai[k - 1] * deg;
    }

/*      ********************************************************* */
/*      * Transformation of the series from ESAPHO for position * */
/*      ********************************************************* */
    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ep4_1.nb[i__ - 1];
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {
      {
        char *str = x + (j + i__ * 27 - 28) * 57;
        char substr[20];
        int k;
        
        for (k = 0; k < 6; k++)
          SCAN(3, d, iar + k);
        
        SCAN(13, lf, &r1);
        SCAN(8, d, ic1);
        for (k = 1; k < 3; k++)
        {
          SCAN(9, d, ic1 + k);
        }
        
        str = y + (j + i__ * 27 - 28) * 51;
        
        SCAN(14, lf, &r2);
        for (k = 0; k < 4; k++)
          SCAN(6, d, ic2 + k);
        SCAN(13, d, ic1 + 3);
      }
      
	    aw = 0.;
	    fw = 0.;
	    for (k = 1; k <= 6; ++k) {
		aw += iar[k - 1] * a0[k - 1];
		fw += iar[k - 1] * f[k - 1];
	    }
	    ep5_1.arg[j + i__ * 31 - 32] = aw;
	    ep5_1.freq[j + i__ * 31 - 32] = fw;
	    ep5_1.freq2[j + i__ * 31 - 32] = (doublereal) (iar[2] + iar[3] + 
		    iar[4]);
	    r1 = r1 + ic1[0] * dnun + ic1[1] * ep1_2.dgam + ic1[2] * ep1_2.de 
		    + ic1[3] * dpn;
	    r2 = r2 + ic2[0] * dnun + ic2[1] * ep1_2.dgam + ic2[2] * ep1_2.de 
		    + ic2[3] * dpn;
	    if (i__ == 1) {
		ep5_1.cs[j + i__ * 31 - 32] = r2 * .001;
		ep5_1.cc[j + i__ * 31 - 32] = r1 * .001;
	    } else {
		ep5_1.cs[j + i__ * 31 - 32] = r1 * .001;
		ep5_1.cc[j + i__ * 31 - 32] = r2 * .001;
	    }
	}
	if (i__ != 3) {
	    for (j = 1; j <= 4; ++j) {
		jp = jmax + j;
		ep5_1.arg[jp + i__ * 31 - 32] = asup[j - 1] * deg;
		ep5_1.freq[jp + i__ * 31 - 32] = fsup[j - 1] * deg;
		ep5_1.freq2[jp + i__ * 31 - 32] = 0.;
		if (i__ == 1) {
		    ep5_1.cs[jp + i__ * 31 - 32] = 0.;
		    ep5_1.cc[jp + i__ * 31 - 32] = rx1s[j - 1] * .001;
		} else {
		    ep5_1.cs[jp + i__ * 31 - 32] = rx2s[j - 1] * .001;
		    ep5_1.cc[jp + i__ * 31 - 32] = 0.;
		}
	    }
	}
    }


/*      ********************************************************* */
/*      * Transformation of the series from ESAPHO for velocity * */
/*      ********************************************************* */
    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ep4_1.nbv[i__ - 1];
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {
      {
        char *str = v + (j + (i__ << 5) - 33) * 61;
        char substr[20];
        int k;
        
        for (k = 0; k < 6; k++)
          SCAN(3, d, iar + k);
        
        SCAN(14, lf, &r1);
        SCAN(10, d, ic1);
        SCAN(9, d, ic1 + 1);
        SCAN(10, d, ic1 + 2);

        str = w + (j + (i__ << 5) - 33) * 56;
        SCAN(15, lf, &r2);
        SCAN(6, d, ic2);
        SCAN(8, d, ic2 + 1);
        SCAN(7, d, ic2 + 2);
        SCAN(6, d, ic2 + 3);
        SCAN(14, d, ic1 + 3);
      }      
	    aw = 0.;
	    fw = 0.;
	    for (k = 1; k <= 6; ++k) {
		aw += iar[k - 1] * a0[k - 1];
		fw += iar[k - 1] * f[k - 1];
	    }
	    ep5_1.argv[j + i__ * 36 - 37] = aw;
	    ep5_1.freqv[j + i__ * 36 - 37] = fw;
	    ep5_1.frev2[j + i__ * 36 - 37] = (doublereal) (iar[2] + iar[3] + 
		    iar[4]);
	    r1 = r1 + ic1[0] * dnun + ic1[1] * ep1_2.dgam + ic1[2] * ep1_2.de 
		    + ic1[3] * dpn;
	    r2 = r2 + ic2[0] * dnun + ic2[1] * ep1_2.dgam + ic2[2] * ep1_2.de 
		    + ic2[3] * dpn;
	    if (i__ != 1) {
		ep5_1.cvs[j + i__ * 36 - 37] = r2 * .001;
		ep5_1.cvc[j + i__ * 36 - 37] = r1 * .001;
	    } else {
		ep5_1.cvs[j + i__ * 36 - 37] = r1 * .001;
		ep5_1.cvc[j + i__ * 36 - 37] = r2 * .001;
	    }
	}
	if (i__ != 3) {
	    for (j = 1; j <= 4; ++j) {
		jp = jmax + j;
		ep5_1.argv[jp + i__ * 36 - 37] = asup[j - 1] * deg;
		ep5_1.freqv[jp + i__ * 36 - 37] = fsup[j - 1] * deg;
		ep5_1.frev2[jp + i__ * 36 - 37] = 0.;
		if (i__ == 1) {
		    ep5_1.cvs[jp + i__ * 36 - 37] = -rx1s[j - 1] * .001 * 
			    ep5_1.argv[jp + i__ * 36 - 37];
		    ep5_1.cvc[j + i__ * 36 - 37] = 0.;
		} else {
		    ep5_1.cvs[jp + i__ * 36 - 37] = 0.;
		    ep5_1.cvc[jp + i__ * 36 - 37] = rx2s[j - 1] * .001 * 
			    ep5_1.argv[jp + i__ * 36 - 37];
		}
	    }
	}
    }


    return 0;
} /* lec_ */

#undef w2p
#undef w1p
#undef y2p
#undef y1p
#undef x2p
#undef x1p
#undef v2p
#undef v1p
#undef w3
#undef w2
#undef w1
#undef v3
#undef y3
#undef y2
#undef y1
#undef x3
#undef x2
#undef x1
#undef v2
#undef v1
#undef y
#undef x
#undef w
#undef v





/* Subroutine */ int lecd_(void)
{
    /* Initialized data */

    static doublereal fipsi = 350.8919885;
    static doublereal fipi = 1.772311e-5;
    static doublereal fid = 284.6378363;
    static doublereal fif = 285.179876;
    static doublereal fil = 285.143868;
    static doublereal filp = .5240207;
    static doublereal fnu = 285.161908;
    static struct {
	char e_1[2183];
	char fill_2[59];
	char e_3[3717];
	char fill_4[767];
	} equiv_9 = { "  0  1  1  0  0  1 23451762.1311-15634508  -727543   "
		" -9779  0  2  2 -1  0  2   -55151.4362    36767        0    "
		"    0  0  0  0  1  0  0    55133.9035   -36755        0     "
		"   0  0  1  1  0 -1  1   -14544.5209     9696      455-35177"
		"119  0  1  1 -2  0  1     5602.8085    -3735   727543       "
		"-2  0  1  1  0  1  1     4890.3869    -3260     -151 1172570"
		"1  0  1  1  0  0  0     4551.1609    -3034        0        0"
		"  0  1  1  0  0  2    -4545.4296     3030        0        0 "
		" 2 -3 -3  0  0 -3    -3251.3962     2167        0        0  "
		"0  3  1  0  0  3    -2659.0343     1772        0        0  0"
		"  1 -1  0  0  1     2592.0639    -1728        0        0  2 "
		"-1 -1  0  0 -1     1868.9622    -1245        0        0  0  "
		"1  1 -1  0  1        0.0000        0        0        0  0  3"
		"  1  0  0  4     -581.0963      387        0        0  0  1 "
		"-1  0  0  2      558.8118     -372        0        0  0  2  "
		"2  0  0  2        0.0000        0        0        0  0  2  0"
		"  1  0  2     -390.0656      260        0        0  0  1  1 "
		" 1  0  1        0.0000        0        0        0  3 -4 -4  "
		"0  0 -4      280.7530     -187        0        0  0  1  1  0"
		"  0 -1      318.5942     -212        0        0  0  1  1  0 "
		" 0  3     -316.7712      211        0        0  0  0  2 -1  "
		"0  0      288.9456     -192        0        0  0  3  1  0  0"
		"  2      255.1685     -170        0        0  0  1 -1  0  0 "
		" 0     -252.9696      168        0        0  1 -2 -2  0  0 -"
		"2      157.2017     -104        0        0  0  3  3 -2  0  3"
		"      148.1044      -98        0        0  3 -2 -2  0  0 -2 "
		"    -126.2662       84        0        0  1  0  0  0  0  0  "
		"   -119.1380       79        0        0  0  0  0  1  0  1   "
		"   125.3480      -83        0        0  0  0  0  1  0 -1    "
		" -118.4993       78        0        0  0  3  1  0  0  5     "
		" -98.7740       65        0        0  0  2  2 -1  0  1      "
		"-94.0167       62        0        0  0  1 -1  0  0  3       "
		"93.6378      -62        0        0  0  2  2 -1  0  3       8"
		"6.4433      -57        0        0  0  2  0 -1  0  2       85"
		".1187      -56        0        0  0  2  0  1  0  3      -84."
		"9226       56        0        0  0  0  2 -1  0 -1       62.9"
		"106      -41        0        0", {0}, "  0  1  1  0  0  1 23"
		"448856.2674-15632570  -727543    -9779  0  0  0  1  0  0    "
		"66507.3861   -44338        0        0  0  2  2 -1  0  2   -5"
		"5124.5490    36749        0        0  0  1  1  0 -1  1   -14"
		"794.1622     9862      455-35177119  0  1  1 -2  0  1     56"
		"52.0343    -3768   727543       -2  0  1  1  0  1  1     488"
		"9.7715    -3259     -151 11725701  0  1  1  0  0  0     4526"
		".5150    -3017        0        0  0  1  1  0  0  2    -4514."
		"4268     3009        0        0  2 -3 -3  0  0 -3     3250.9"
		"885    -2167        0        0  0  1 -1  0  0  1    -2843.19"
		"94     1895        0        0  0  3  1  0  0  3    -2653.737"
		"5     1769        0        0  2 -1 -1  0  0 -1     1868.6828"
		"    -1245        0        0  0  0  0  0  0  0        0.0000 "
		"       0        0        0  0  1  1 -1  0  1        0.0000  "
		"      0        0        0  0  1 -1  0  0  2     -632.9468   "
		"   421        0        0  0  3  1  0  0  4     -579.9308    "
		"  386        0        0  0  2  0  1  0  2     -399.2736     "
		" 266        0        0  0  2  2  0  0  2        0.0000      "
		"  0        0        0  0  1  1  1  0  1        0.0000       "
		" 0        0        0  3 -4 -4  0  0 -4     -280.7182      18"
		"7        0        0  0  1  1  0  0 -1      316.9349     -211"
		"        0        0  0  1  1  0  0  3     -314.5569      209 "
		"       0        0  0  0  2 -1  0  0      294.7831     -196  "
		"      0        0  0  1 -1  0  0  0      269.2981     -179   "
		"     0        0  0  3  1  0  0  2      254.6917     -169    "
		"    0        0  1 -2 -2  0  0 -2     -157.1732      104     "
		"   0        0  0  3  3 -2  0  3      148.0042      -98      "
		"  0        0  3 -2 -2  0  0 -2     -126.2468       84       "
		" 0        0  1  0  0  0  0  0     -119.1139       79        "
		"0        0  0  1 -1  0  0  3     -109.5100       73        0"
		"        0  0  2  0 -1  0  2      106.0245      -70        0 "
		"       0  0  3  1  0  0  5      -98.5724       65        0  "
		"      0  0  2  2 -1  0  1      -93.6819       62        0   "
		"     0  0  2  0  1  0  3      -86.9127       57        0    "
		"    0  0  2  2 -1  0  3       86.0623      -57        0     "
		"   0  0  0  0  1  0  1       70.1733      -46        0      "
		"  0  0  0  0  1  0 -1      -66.2919       44        0       "
		" 0  0  0  2 -1  0 -1       64.1928      -42        0        "
		"0  0  0  0  1  0  0   726494.5032  -484329 46897178     -303"
		"  0  1  1  0  0  1  -367660.3813   245106        0        0 "
		" 0  1 -1  0  0  1    -6081.1262     4054        0        0  "
		"0  2  2 -1  0  2     2149.8425    -1433        0        0  0"
		"  1  1  0  0  2     1844.8540    -1229        0        0  0 "
		" 1  1 -2  0  1    -1799.3963     1199        0        0  0  "
		"1  1  0  0  0    -1729.1650     1152        0        0  0  1"
		" -1  0  0  2    -1339.1617      892        0        0  0  1 "
		"-1  0  0  0      574.0342     -382        0        0  0  0  "
		"0  1 -1  0     -454.8726      303   -29334 -1091184  0  0  2"
		" -1  0  0      417.5738     -278        0        0  0  3  1 "
		" 0  0  3      333.6614     -222        0        0  0  1  1  "
		"0 -1  1      230.5734     -153        0        0  0  1 -1  0"
		"  0  3     -228.7458      152        0        0  0  0  0  1 "
		" 1  0      151.3386     -100     9778   363728  0  1  1  0  "
		"0  3      126.3330      -84        0        0  0  1  1  0  0"
		" -1     -123.0024       82        0        0  2 -2 -2 -1  0 "
		"-2      100.4703      -66        0        0  0  0  2 -1  0 -"
		"1       92.7574      -61        0        0  0  1  1  0  1  1"
		"      -77.0433       51        0        0  0  3  1  0  0  4 "
		"      72.3397      -48        0        0  2 -2 -2  1  0 -2  "
		"     63.1447      -42        0        0  0  0  0  0  0  0   "
		"     0.0000        0        0        0  0  1  1 -1  0  1    "
		"    0.0000        0        0        0  2 -3 -3  0  0 -3     "
		" -51.4568       34        0        0" };

    static struct {
	char e_1[370];
	char fill_2[10];
	char e_3[630];
	char fill_4[130];
	} equiv_11 = { "    0.0000    0.0000    0.0000    0.0000    0.0000  "
		"  0.0000    0.0000    0.0000   -0.0020    0.0000    0.0000  "
		"  0.0007-1157.0034    0.0000    0.0000  395.4635    0.0000 -"
		"383.9143   45.6488    0.0000    0.0000    0.0000    0.0000  "
		"  0.0000  -11.6027    0.0000  -20.5254    8.8034    0.0000  "
		"  0.0000    0.0000    0.0000    0.0000    0.0000    0.0000  "
		"  0.0000    0.0000", {0}, "    0.0000    0.0000    0.0000   "
		" 0.0000    0.0000    0.0000    0.0000    0.0000   -0.0020   "
		" 0.0000    0.0000   -0.0007-1180.2676 1156.5679    0.0000   "
		" 0.0000    0.0000 -395.4136  383.6747   45.6431    0.0000   "
		" 0.0000    0.0000    0.0000    0.0000  -11.6006    0.0000   "
		"20.5223   -8.8016    0.0000    0.0000    0.0000    0.0000   "
		" 0.0000    0.0000    0.0000    0.0000    0.0000    0.0000   "
		" 0.0000    0.0000    0.0000    0.0000    0.0000    0.0000   "
		" 0.0000    0.0000    0.0000    0.0000    0.0000    0.0000   "
		" 0.0000    0.0000    0.0000    0.0000    0.0000    0.0000   "
		" 0.0000    0.0000    0.0000   57.4788  -54.7313    0.0000" };

    static struct {
	char e_1[4410];
	char fill_2[693];
	} equiv_20 = { "  0  1  1  0  0  1 -116719753.045 -38906584   362099"
		"1     48670  0  2  2 -1  0  2     274472.052     91511      "
		"   0         0  0  0  0  1  0  0    -274419.644    -91460   "
		"      0         0  0  1  1 -2  0  1      27888.799      9293"
		"   3621451        -9  0  1  1  0  1  1     -48677.506    -16"
		"229      1503-116714258  0  1  1  0  0  0     -22609.568    "
		" -7578         0         0  0  1  1  0  0  2      22664.240 "
		"     7514         0         0  2 -3 -3  0  0 -3      -8722.2"
		"00    -42733         0         0  0  3  1  0  0  3      1328"
		"2.692      4382         0         0  0  1 -1  0  0  1      1"
		"2853.324      4332         0         0  2 -1 -1  0  0 -1    "
		" -13590.022     18354         0         0  0  3  1  0  0  4 "
		"      2908.069       955         0         0  0  1 -1  0  0 "
		" 2       2765.881       939         0         0  0  2  2  0 "
		" 0  2          0.000         0         0         0  0  2  0 "
		" 1  0  2       1948.621       642         0         0  0  1 "
		" 1  1  0  1          0.000         0         0         0  0 "
		" 1  1  0  0 -1      -1579.820      -534         0         0 "
		" 0  1  1  0  0  3       1582.369       522         0        "
		" 0  0  0  2 -1  0  0      -1432.710      -486         0     "
		"    0  0  3  1  0  0  2      -1272.311      -422         0  "
		"       0  0  1 -1  0  0  0      -1256.720      -424         "
		"0         0  1 -2 -2  0  0 -2        602.052      1166      "
		"   0         0  0  3  3 -2  0  3       -737.024      -249   "
		"      0         0  3 -2 -2  0  0 -2       1062.989     -1964"
		"         0         0  1  0  0  0  0  0        729.627      -"
		"483         0         0  0  0  0  1  0  1       -625.045    "
		"  -209         0         0  0  0  0  1  0 -1        588.726 "
		"      202         0         0  0  1  1  0  0  1  116705290.5"
		"20  38901767  -3620991    -48670  0  0  0  1  0  0     33102"
		"9.222    110323         0         0  0  2  2 -1  0  2    -27"
		"4338.243    -91466         0         0  0  1  1 -2  0  1    "
		" -28133.827     -9374  -3621451         9  0  1  1  0  1  1 "
		"     48671.381     16233     -1503 116714258  0  1  1  0  0 "
		" 0      22487.130      7540         0         0  0  1  1  0 "
		" 0  2     -22509.655     -7465         0         0  2 -3 -3 "
		" 0  0 -3      -8721.107    -42727         0         0  0  1 "
		"-1  0  0  1      14098.635      4753         0         0  0 "
		" 3  1  0  0  3     -13256.232     -4370         0         0 "
		" 2 -1 -1  0  0 -1      13587.990    -18353         0        "
		" 0  0  1 -1  0  0  2       3132.818      1066         0     "
		"    0  0  3  1  0  0  4      -2902.236      -954         0  "
		"       0  0  2  0  1  0  2      -1994.620      -658         "
		"0         0  0  2  2  0  0  2          0.000         0      "
		"   0         0  0  1  1  1  0  1          0.000         0   "
		"      0         0  0  1  1  0  0 -1       1571.592       531"
		"         0         0  0  1  1  0  0  3      -1571.308      -"
		"521         0         0  0  0  2 -1  0  0       1461.654    "
		"   495         0         0  0  1 -1  0  0  0      -1337.837 "
		"     -451         0         0  0  3  1  0  0  2       1269.9"
		"34       424         0         0  1 -2 -2  0  0 -2        60"
		"1.943      1166         0         0  0  3  3 -2  0  3       "
		" 736.525       248         0         0  3 -2 -2  0  0 -2    "
		"  -1062.825      1963         0         0  1  0  0  0  0  0 "
		"      -729.480       483         0         0  0  1 -1  0  0 "
		" 3        541.026       184         0         0  0  2  0 -1 "
		" 0  2       -525.780      -180         0         0  0  0  0 "
		" 1  0  0    3616003.041   1205108 233422741     -1508  0  1 "
		" 1  0  0  1   -1829850.937   -609954         0         0  0 "
		" 1 -1  0  0  1      30154.613     10163         0         0 "
		" 0  2  2 -1  0  2      10699.117      3568         0        "
		" 0  0  1  1  0  0  2       9198.737      3053         0     "
		"    0  0  1  1 -2  0  1       8956.758      2987         0  "
		"       0  0  1  1  0  0  0      -8590.264     -2883         "
		"0         0  0  1 -1  0  0  2       6628.283      2250      "
		"   0         0  0  1 -1  0  0  0      -2851.726      -959   "
		"      0         0  0  0  2 -1  0  0       2070.500       699"
		"         0         0  0  3  1  0  0  3       1666.741       "
		"551         0         0  0  1 -1  0  0  3       1130.103    "
		"   387         0         0  0  0  0  1  1  0       1506.429 "
		"      511     97330   3620558  0  1  1  0  0  3        631.0"
		"72       209         0         0  0  1  1  0  0 -1       -60"
		"9.935      -205         0         0  0  1  1  0  1  1       "
		"-766.867      -259         0         0" };

    static struct {
	char e_1[770];
	char fill_2[121];
	} equiv_21 = { "      0.000      0.000      0.000      0.000      0."
		"000      0.000      0.000      0.005      0.000      0.000  "
		"    0.005      0.000      0.000   3936.455      0.000  -3821"
		".615      0.000      0.000      0.000      0.000      0.000 "
		"    44.436      0.000   -172.796     53.914      0.000      "
		"0.000      0.000      0.000      0.000      0.000      0.000"
		"      0.000      0.000     -0.005      0.000      0.000     "
		" 0.005      0.000      0.000      0.000   3935.958  -3819.23"
		"0      0.000      0.000      0.000      0.000      0.000    "
		"-44.428      0.000   -172.770     53.903      0.000      0.0"
		"00      0.000      0.000      0.000      0.000      0.000   "
		"   0.000      0.000      0.000      0.000      0.000      0."
		"000      0.000      0.000      0.000      0.000      0.000" };


    /* System generated locals */
    integer i__1;
    icilist ici__1;

    /* Local variables */
    static doublereal f[6];
    static integer i__, j, k;
#define v ((char *)&equiv_20)
#define w ((char *)&equiv_21)
#define x ((char *)&equiv_9)
#define y ((char *)&equiv_11)
    static doublereal a0[6], r1, r2;
#define v1 ((char *)&equiv_20)
#define v2 ((char *)&equiv_20 + 1701)
#define x1 ((char *)&equiv_9)
#define x2 ((char *)&equiv_9 + 2242)
#define x3 ((char *)&equiv_9 + 4484)
#define y1 ((char *)&equiv_11)
#define y2 ((char *)&equiv_11 + 380)
#define y3 ((char *)&equiv_11 + 760)
#define v3 ((char *)&equiv_20 + 3402)
#define w1 ((char *)&equiv_21)
#define w2 ((char *)&equiv_21 + 297)
#define w3 ((char *)&equiv_21 + 594)
    static doublereal ai[6], pf[6], aw, fw;
    static integer ic1[3];
#define v1p ((char *)&equiv_20 + 1134)
#define v2p ((char *)&equiv_20 + 2835)
#define x1p ((char *)&equiv_9 + 1062)
#define x2p ((char *)&equiv_9 + 3363)
#define x3p ((char *)&equiv_9 + 5546)
#define y1p ((char *)&equiv_11 + 180)
#define y2p ((char *)&equiv_11 + 570)
#define y3p ((char *)&equiv_11 + 940)
#define w1p ((char *)&equiv_21 + 198)
#define w2p ((char *)&equiv_21 + 495)
    static doublereal deg;
    static integer iar[6], jmax;
    static doublereal dnun;

/*      ********************************************************* */
/*      * Transformation of the series from ESADE.              * */
/*      * Introduction of the new constants in the coefficients * */
/*      * and in the mean motions                               * */
/*      ********************************************************* */







/*      *************************** */
/*      * Coefficients from ESADE * */
/*      *************************** */

/*      ************************ */
/*      * Arguments in J2000.0 * */
/*      ************************ */
    ai[0] = ep3_1.psi0;
    ai[1] = ep3_1.pipet;
    ai[2] = ed1_2.al[0] - ep3_1.pipet - ep3_1.alm;
    ai[3] = ed1_2.al[0] - ed1_2.al[1];
    ai[4] = ed1_2.al[0] - ed1_2.al[2];
    ai[5] = ep3_1.alm;

/*      *************************** */
/*      * Mean motions from ESADE * */
/*      *************************** */

    dnun = ed1_2.dnu / fnu;
    pf[0] = fipsi;
    pf[1] = fipi;
    pf[2] = fid;
    pf[3] = fif;
    pf[4] = fil;
    pf[5] = filp;

/*      ************************* */
/*      * Conversion to radians * */
/*      ************************* */
    deg = .017453292519943292;
    for (k = 1; k <= 6; ++k) {
	f[k - 1] = pf[k - 1] * deg;
	a0[k - 1] = ai[k - 1] * deg;
    }

/*      ******************************************************** */
/*      * Transformation of the series from ESADE for position * */
/*      ******************************************************** */
    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ed4_2.nb[i__ - 1];
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {
    
      {
        char *str = x + (j + i__ * 38 - 39) * 59;
        char substr[20];
        int k;
        
        for (k = 0; k < 6; k++)
          SCAN(3, d, iar + k);
        
        SCAN(14, lf, &r1);
        for (k = 0; k < 3; k++)
          SCAN(9, d, ic1 + k);
        
        str = y + (j + i__ * 38 - 39) * 10;
        SCAN(10, lf, &r2);
      }
      
      
	    aw = 0.;
	    fw = 0.;
	    for (k = 1; k <= 6; ++k) {
		aw += iar[k - 1] * a0[k - 1];
		fw += iar[k - 1] * f[k - 1];
	    }
	    ed5_1.arg[j + i__ * 38 - 39] = aw;
	    ed5_1.freq[j + i__ * 38 - 39] = fw;
	    ed5_1.freq2[j + i__ * 38 - 39] = (doublereal) (iar[2] + iar[3] + 
		    iar[4]);
	    r1 = r1 + ic1[0] * dnun + ic1[1] * ed1_2.dgam + ic1[2] * ed1_2.de;
	    r2 = r2;
	    if (i__ == 1) {
		ed5_1.cs[j + i__ * 38 - 39] = r2 * .001;
		ed5_1.cc[j + i__ * 38 - 39] = r1 * .001;
	    } else {
		ed5_1.cs[j + i__ * 38 - 39] = r1 * .001;
		ed5_1.cc[j + i__ * 38 - 39] = r2 * .001;
	    }
	}
    }


/*      ********************************************************* */
/*      * Transformation of the series from ESADE for velocity  * */
/*      ********************************************************* */
    for (i__ = 1; i__ <= 3; ++i__) {
	jmax = ed4_2.nbv[i__ - 1];
	i__1 = jmax;
	for (j = 1; j <= i__1; ++j) {

      {
        char *str = v + (j + i__ * 27 - 28) * 63;
        char substr[20];
      
        for (k = 0; k < 6; k++)
          SCAN(3, d, iar + k);
        
        SCAN(15, lf, &r1);
        for (k = 0; k < 3; k++)
          SCAN(10, d, ic1 + k);
        
        str = w + (j + i__ * 27 - 28) * 11;
        SCAN(11, lf, &r2);
      }    
      
	    aw = 0.;
	    fw = 0.;
	    for (k = 1; k <= 6; ++k) {
		aw += iar[k - 1] * a0[k - 1];
		fw += iar[k - 1] * f[k - 1];
	    }
	    ed5_1.argv[j + i__ * 27 - 28] = aw;
	    ed5_1.freqv[j + i__ * 27 - 28] = fw;
	    ed5_1.frev2[j + i__ * 27 - 28] = (doublereal) (iar[2] + iar[3] + 
		    iar[4]);
	    r1 = r1 + ic1[0] * dnun + ic1[1] * ed1_2.dgam + ic1[2] * ed1_2.de;
	    r2 = r2;
	    if (i__ != 1) {
		ed5_1.cvs[j + i__ * 27 - 28] = r2 * .001;
		ed5_1.cvc[j + i__ * 27 - 28] = r1 * .001;
	    } else {
		ed5_1.cvs[j + i__ * 27 - 28] = r1 * .001;
		ed5_1.cvc[j + i__ * 27 - 28] = r2 * .001;
	    }
	}
    }


    return 0;
} /* lecd_ */

#undef w2p
#undef w1p
#undef y3p
#undef y2p
#undef y1p
#undef x3p
#undef x2p
#undef x1p
#undef v2p
#undef v1p
#undef w3
#undef w2
#undef w1
#undef v3
#undef y3
#undef y2
#undef y1
#undef x3
#undef x2
#undef x1
#undef v2
#undef v1
#undef y
#undef x
#undef w
#undef v
