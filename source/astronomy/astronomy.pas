unit astronomy;
(*
    Astronomy Engine for C/C++.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2023 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*)


{$MINENUMSIZE 4}

interface

const
  {$IF Defined(WIN32)}
  LIB_ASTRONOMY = 'astronomy.dll';
  _PU = '';
  {$ELSE}
    {$MESSAGE Error 'Unsupported platform'}
  {$ENDIF}

const
  // The speed of light in AU/day.
  C_AUDAY = 173.1446326846693;
  // The number of kilometers in one astronomical unit (AU).
  KM_PER_AU = 1.4959787069098932e+8;
  // The number of astronomical units per light-year.
  AU_PER_LY = 63241.07708807546;
  // The factor to convert degrees to radians = pi/180.
  DEG2RAD = 0.017453292519943296;
  // The factor to convert sidereal hours to radians = pi/12.
  HOUR2RAD = 0.2617993877991494365;
  // The factor to convert radians to degrees = 180/pi.
  RAD2DEG = 57.295779513082321;
  // The factor to convert radians to sidereal hours = 12/pi.
  RAD2HOUR = 3.819718634205488;
  // The mean radius of the Sun's photosphere, expressed in kilometers. The Sun is nearly spherical.
  SUN_RADIUS_KM = 695700.0;
  // The equatorial radius of Mercury, expressed in kilometers.
  MERCURY_EQUATORIAL_RADIUS_KM = 2440.5;
  // The polar radius of Mercury, expressed in kilometers.
  MERCURY_POLAR_RADIUS_KM = 2438.3;
  // The mean radius of Venus, expressed in kilometers. Venus is nearly spherical.
  VENUS_RADIUS_KM = 6051.8;
  // The equatorial radius of the Earth, expressed in kilometers.
  EARTH_EQUATORIAL_RADIUS_KM = 6378.1366;
  // The Earth's polar radius divided by its equatorial radius.
  EARTH_FLATTENING = 0.996647180302104;
  // The polar radius of the Earth, expressed in kilometers.
  EARTH_POLAR_RADIUS_KM = (EARTH_EQUATORIAL_RADIUS_KM*EARTH_FLATTENING);
  // The equatorial radius of the Moon, expressed in kilometers.
  MOON_EQUATORIAL_RADIUS_KM = 1738.1;
  // The polar radius of the Moon, expressed in kilometers.
  MOON_POLAR_RADIUS_KM = 1736.0;
  // The equatorial radius of Mars, expressed in kilometers.
  MARS_EQUATORIAL_RADIUS_KM = 3396.2;
  // The polar radius of Mars, expressed in kilometers.
  MARS_POLAR_RADIUS_KM = 3376.2;
(*
    Jupiter radius data are nominal values obtained from:
    https://www.iau.org/static/resolutions/IAU2015_English.pdf
    https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
*)
  // The equatorial radius of Jupiter, expressed in kilometers.
  JUPITER_EQUATORIAL_RADIUS_KM = 71492.0;
  // The polar radius of Jupiter, expressed in kilometers.
  JUPITER_POLAR_RADIUS_KM = 66854.0;
  // The volumetric mean radius of Jupiter, expressed in kilometers.
  JUPITER_MEAN_RADIUS_KM = 69911.0;
(*
    The radii of Jupiter's four major moons are obtained from:
    https://ssd.jpl.nasa.gov/?sat_phys_par
*)
  // The mean radius of Jupiter's moon Io, expressed in kilometers.
  IO_RADIUS_KM = 1821.6;
  // The mean radius of Jupiter's moon Europa, expressed in kilometers.
  EUROPA_RADIUS_KM = 1560.8;
  // The mean radius of Jupiter's moon Ganymede, expressed in kilometers.
  GANYMEDE_RADIUS_KM = 2631.2;
  // The mean radius of Jupiter's moon Callisto, expressed in kilometers.
  CALLISTO_RADIUS_KM = 2410.3;
  // The equatorial radius of Saturn, expressed in kilometers.
  SATURN_EQUATORIAL_RADIUS_KM = 60268.0;
  // The polar radius of Saturn, expressed in kilometers.
  SATURN_POLAR_RADIUS_KM = 54364.0;
  // The equatorial radius of Uranus, expressed in kilometers.
  URANUS_EQUATORIAL_RADIUS_KM = 25559.0;
  // The polar radius of Uranus, expressed in kilometers.
  URANUS_POLAR_RADIUS_KM = 24973.0;
  // The equatorial radius of Neptune, expressed in kilometers.
  NEPTUNE_EQUATORIAL_RADIUS_KM = 24764.0;
  // The polar radius of Neptune, expressed in kilometers.
  NEPTUNE_POLAR_RADIUS_KM = 24341.0;
  // The mean radius of Pluto, expressed in kilometers. Pluto is nearly spherical.
  PLUTO_RADIUS_KM = 1188.3;
  TIME_TEXT_BYTES = 28;
  { TODO : Unable to convert function-like macro: }
  (* Astronomy_LongitudeFromSun ( body , time ) ( Astronomy_PairLongitude ( ( body ) , BODY_SUN , ( time ) ) ) *)
  { TODO : Unable to convert function-like macro: }
  (* Astronomy_SearchHourAngle ( body , observer , hourAngle , startTime ) Astronomy_SearchHourAngleEx ( ( body ) , ( observer ) , ( hourAngle ) , ( startTime ) , + 1 ) *)
  { TODO : Unable to convert function-like macro: }
  (* Astronomy_SearchRiseSet ( body , observer , direction , startTime , limitDays ) Astronomy_SearchRiseSetEx ( ( body ) , ( observer ) , ( direction ) , ( startTime ) , ( limitDays ) , 0.0 ) *)

type
  // Forward declarations
  Pastro_grav_sim_s = Pointer;
  PPastro_grav_sim_s = ^Pastro_grav_sim_s;
  Pastro_time_t = ^astro_time_t;
  Pastro_utc_t = ^astro_utc_t;
  Pastro_vector_t = ^astro_vector_t;
  Pastro_state_vector_t = ^astro_state_vector_t;
  Pastro_spherical_t = ^astro_spherical_t;
  Pastro_angle_result_t = ^astro_angle_result_t;
  Pastro_observer_t = ^astro_observer_t;
  Pastro_equatorial_t = ^astro_equatorial_t;
  Pastro_ecliptic_t = ^astro_ecliptic_t;
  Pastro_horizon_t = ^astro_horizon_t;
  Pastro_rotation_t = ^astro_rotation_t;
  Pastro_atmosphere_t = ^astro_atmosphere_t;
  Pastro_search_result_t = ^astro_search_result_t;
  Pastro_seasons_t = ^astro_seasons_t;
  Pastro_moon_quarter_t = ^astro_moon_quarter_t;
  Pastro_func_result_t = ^astro_func_result_t;
  Pastro_elongation_t = ^astro_elongation_t;
  Pastro_hour_angle_t = ^astro_hour_angle_t;
  Pastro_illum_t = ^astro_illum_t;
  Pastro_apsis_t = ^astro_apsis_t;
  Pastro_lunar_eclipse_t = ^astro_lunar_eclipse_t;
  Pastro_global_solar_eclipse_t = ^astro_global_solar_eclipse_t;
  Pastro_eclipse_event_t = ^astro_eclipse_event_t;
  Pastro_local_solar_eclipse_t = ^astro_local_solar_eclipse_t;
  Pastro_transit_t = ^astro_transit_t;
  Pastro_constellation_t = ^astro_constellation_t;
  Pastro_libration_t = ^astro_libration_t;
  Pastro_axis_t = ^astro_axis_t;
  Pastro_jupiter_moons_t = ^astro_jupiter_moons_t;
  Pastro_node_event_t = ^astro_node_event_t;

  (**
   * @brief Indicates success/failure of an Astronomy Engine function call.
   *)
  astro_status_t = (
    (** The operation was successful. *)
    ASTRO_SUCCESS = 0,
    (** A placeholder that can be used for data that is not yet initialized. *)
    ASTRO_NOT_INITIALIZED = 1,
    (** The celestial body was not valid. Different sets of bodies are supported depending on the function. *)
    ASTRO_INVALID_BODY = 2,
    (** A numeric solver failed to converge. This should not happen unless there is a bug in Astronomy Engine. *)
    ASTRO_NO_CONVERGE = 3,
    (** The provided date/time is outside the range allowed by this function. *)
    ASTRO_BAD_TIME = 4,
    (** Vector magnitude is too small to be normalized into a unit vector. *)
    ASTRO_BAD_VECTOR = 5,
    (** Search was not able to find an ascending root crossing of the function in the specified time interval. *)
    ASTRO_SEARCH_FAILURE = 6,
    (** The Earth cannot be treated as a celestial body seen from an observer on the Earth itself. *)
    ASTRO_EARTH_NOT_ALLOWED = 7,
    (** No lunar quarter occurs inside the specified time range. *)
    ASTRO_NO_MOON_QUARTER = 8,
    (** Internal error: Astronomy_NextMoonQuarter found the wrong moon quarter. *)
    ASTRO_WRONG_MOON_QUARTER = 9,
    (** A self-check failed inside the code somewhere, indicating a bug needs to be fixed. *)
    ASTRO_INTERNAL_ERROR = 10,
    (** A parameter value passed to a function was not valid. *)
    ASTRO_INVALID_PARAMETER = 11,
    (** Special-case logic for finding Neptune/Pluto apsis failed. *)
    ASTRO_FAIL_APSIS = 12,
    (** A provided buffer's size is too small to receive the requested data. *)
    ASTRO_BUFFER_TOO_SMALL = 13,
    (** An attempt to allocate memory failed. *)
    ASTRO_OUT_OF_MEMORY = 14,
    (** The provided initial state vectors did not have matching times. *)
    ASTRO_INCONSISTENT_TIMES = 15);
  Pastro_status_t = ^astro_status_t;

  (**
   * @brief A date and time used for astronomical calculations.
   *
   * This type is of fundamental importance to Astronomy Engine.
   * It is used to represent dates and times for all astronomical calculations.
   * It is also included in the values returned by many Astronomy Engine functions.
   *
   * To create a valid astro_time_t value from scratch, call #Astronomy_MakeTime
   * (for a given calendar date and time) or #Astronomy_CurrentTime (for the system's
   * current date and time).
   *
   * To adjust an existing astro_time_t by a certain real number of days,
   * call #Astronomy_AddDays.
   *
   * The astro_time_t type contains `ut` to represent Universal Time (UT1/UTC) and
   * `tt` to represent Terrestrial Time (TT, also known as *ephemeris time* ).
   * The difference `tt-ut` is known as *&Delta;T*, using a best-fit piecewise model devised by
   * [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).
   *
   * Both `tt` and `ut` are necessary for performing different astronomical calculations.
   * Indeed, certain calculations (such as rise/set times) require both time scales.
   * See the documentation for the `ut` and `tt` fields for more detailed information.
   *
   * In cases where `astro_time_t` is included in a structure returned by
   * a function that can fail, the `astro_status_t` field `status` will contain a value
   * other than `ASTRO_SUCCESS`; in that case the `ut` and `tt` will hold `NAN` (not a number).
   * In general, when there is an error code stored in a struct field `status`, the
   * caller should ignore all other values in that structure, including the `ut` and `tt`
   * inside `astro_time_t`.
   *)
  astro_time_t = record
    (**
     * @brief   UT1/UTC number of days since noon on January 1, 2000.
     *
     * The floating point number of days of Universal Time since noon UTC January 1, 2000.
     * Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
     * not exactly equivalent; UTC and UT1 can disagree by up to &plusmn;0.9 seconds.
     * This approximation is sufficient for the accuracy requirements of Astronomy Engine.
     *
     * Universal Time Coordinate (UTC) is the international standard for legal and civil
     * timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
     * UTC is kept in sync with unpredictable observed changes in the Earth's rotation
     * by occasionally adding leap seconds as needed.
     *
     * UT1 is an idealized time scale based on observed rotation of the Earth, which
     * gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
     * large scale weather events like hurricanes, and internal seismic and convection effects.
     * Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
     * is adjusted by a scheduled whole number of leap seconds as needed.
     *
     * The value in `ut` is appropriate for any calculation involving the Earth's rotation,
     * such as calculating rise/set times, culumination, and anything involving apparent
     * sidereal time.
     *
     * Before the era of atomic timekeeping, days based on the Earth's rotation
     * were often known as *mean solar days*.
     *)
    ut: Double;
    (**
     * @brief   Terrestrial Time days since noon on January 1, 2000.
     *
     * Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
     * In this system, days are not based on Earth rotations, but instead by
     * the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
     * divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
     * for changes in the Earth's rotation.
     *
     * The value in `tt` is used for calculations of movements not involving the Earth's rotation,
     * such as the orbits of planets around the Sun, or the Moon around the Earth.
     *
     * Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
     *)
    tt: Double;
    (**
     * @brief   For internal use only. Used to optimize Earth tilt calculations.
     *)
    psi: Double;
    (**
     * @brief   For internal use only.  Used to optimize Earth tilt calculations.
     *)
    eps: Double;
    (**
     * @brief   For internal use only.  Lazy-caches sidereal time (Earth rotation).
     *)
    st: Double;
  end;

  (**
   * @brief A calendar date and time expressed in UTC.
   *)
  astro_utc_t = record
    (** The year value, e.g. 2019. *)
    year: Integer;
    (** The month value: 1=January, 2=February, ..., 12=December. *)
    month: Integer;
    (** The day of the month in the range 1..31. *)
    day: Integer;
    (** The hour of the day in the range 0..23. *)
    hour: Integer;
    (** The minute of the hour in the range 0..59. *)
    minute: Integer;
    (** The floating point number of seconds in the range [0,60). *)
    second: Double;
  end;

  (**
   * @brief A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).
   *)
  astro_vector_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The Cartesian x-coordinate of the vector in AU. *)
    x: Double;
    (** The Cartesian y-coordinate of the vector in AU. *)
    y: Double;
    (** The Cartesian z-coordinate of the vector in AU. *)
    z: Double;
    (** The date and time at which this vector is valid. *)
    t: astro_time_t;
  end;

  (**
   * @brief A state vector that contains a position (AU) and velocity (AU/day).
   *)
  astro_state_vector_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The Cartesian position x-coordinate of the vector in AU. *)
    x: Double;
    (** The Cartesian position y-coordinate of the vector in AU. *)
    y: Double;
    (** The Cartesian position z-coordinate of the vector in AU. *)
    z: Double;
    (** The Cartesian velocity x-coordinate of the vector in AU/day. *)
    vx: Double;
    (** The Cartesian velocity y-coordinate of the vector in AU/day. *)
    vy: Double;
    (** The Cartesian velocity z-coordinate of the vector in AU/day. *)
    vz: Double;
    (** The date and time at which this state vector is valid. *)
    t: astro_time_t;
  end;

  (**
   * @brief Spherical coordinates: latitude, longitude, distance.
   *)
  astro_spherical_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The latitude angle: -90..+90 degrees. *)
    lat: Double;
    (** The longitude angle: 0..360 degrees. *)
    lon: Double;
    (** Distance in AU. *)
    dist: Double;
  end;

  (**
   * @brief An angular value expressed in degrees.
   *)
  astro_angle_result_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** An angle expressed in degrees. *)
    angle: Double;
  end;

  (**
   * @brief A celestial body.
   *)
  astro_body_t = (
    (** An invalid or undefined celestial body. *)
    BODY_INVALID = -1,
    (** Mercury *)
    BODY_MERCURY = 0,
    (** Venus *)
    BODY_VENUS = 1,
    (** Earth *)
    BODY_EARTH = 2,
    (** Mars *)
    BODY_MARS = 3,
    (** Jupiter *)
    BODY_JUPITER = 4,
    (** Saturn *)
    BODY_SATURN = 5,
    (** Uranus *)
    BODY_URANUS = 6,
    (** Neptune *)
    BODY_NEPTUNE = 7,
    (** Pluto *)
    BODY_PLUTO = 8,
    (** Sun *)
    BODY_SUN = 9,
    (** Moon *)
    BODY_MOON = 10,
    (** Earth/Moon Barycenter *)
    BODY_EMB = 11,
    (** Solar System Barycenter *)
    BODY_SSB = 12,
    (** user-defined star #1 *)
    BODY_STAR1 = 101,
    (** user-defined star #2 *)
    BODY_STAR2 = 102,
    (** user-defined star #3 *)
    BODY_STAR3 = 103,
    (** user-defined star #4 *)
    BODY_STAR4 = 104,
    (** user-defined star #5 *)
    BODY_STAR5 = 105,
    (** user-defined star #6 *)
    BODY_STAR6 = 106,
    (** user-defined star #7 *)
    BODY_STAR7 = 107,
    (** user-defined star #8 *)
    BODY_STAR8 = 108);
  Pastro_body_t = ^astro_body_t;

  (**
   * @brief The location of an observer on (or near) the surface of the Earth.
   *
   * This structure is passed to functions that calculate phenomena as observed
   * from a particular place on the Earth.
   *
   * You can create this structure directly, or you can call the convenience function
   * #Astronomy_MakeObserver to create one for you.
   *)
  astro_observer_t = record
    (** Geographic latitude in degrees north (positive) or south (negative) of the equator. *)
    latitude: Double;
    (** Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England. *)
    longitude: Double;
    (** The height above (positive) or below (negative) sea level, expressed in meters. *)
    height: Double;
  end;

  (**
   * @brief Equatorial angular and cartesian coordinates.
   *
   * Coordinates of a celestial body as seen from the Earth (geocentric or topocentric, depending on context),
   * oriented with respect to the projection of the Earth's equator onto the sky.
   *)
  astro_equatorial_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** right ascension in sidereal hours. *)
    ra: Double;
    (** declination in degrees *)
    dec: Double;
    (** distance to the celestial body in AU. *)
    dist: Double;
    (** equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north. *)
    vec: astro_vector_t;
  end;

  (**
   * @brief Ecliptic angular and Cartesian coordinates.
   *
   * Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
   * oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).
   *)
  astro_ecliptic_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** Cartesian ecliptic vector: x=equinox, y=90 degrees prograde in ecliptic plane, z=northward perpendicular to ecliptic. *)
    vec: astro_vector_t;
    (** Latitude in degrees north (positive) or south (negative) of the ecliptic plane. *)
    elat: Double;
    (** Longitude in degrees around the ecliptic plane prograde from the equinox. *)
    elon: Double;
  end;

  (**
   * @brief Coordinates of a celestial body as seen by a topocentric observer.
   *
   * Contains horizontal and equatorial coordinates seen by an observer on or near
   * the surface of the Earth (a topocentric observer).
   * Optionally corrected for atmospheric refraction.
   *)
  astro_horizon_t = record
    (** Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West. *)
    azimuth: Double;
    (** Angle in degrees above (positive) or below (negative) the observer's horizon. *)
    altitude: Double;
    (** Right ascension in sidereal hours. *)
    ra: Double;
    (** Declination in degrees. *)
    dec: Double;
  end;

  (**
   * @brief Contains a rotation matrix that can be used to transform one coordinate system to another.
   *)
  astro_rotation_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** A normalized 3x3 rotation matrix. *)
    rot: array [0..2] of array [0..2] of Double;
  end;

  (**
   * @brief Selects whether to correct for atmospheric refraction, and if so, how.
   *)
  astro_refraction_t = (
    (** No atmospheric refraction correction (airless). *)
    REFRACTION_NONE = 0,
    (** Recommended correction for standard atmospheric refraction. *)
    REFRACTION_NORMAL = 1,
    (** Used only for compatibility testing with JPL Horizons online tool. *)
    REFRACTION_JPLHOR = 2);
  Pastro_refraction_t = ^astro_refraction_t;

  (**
   * @brief Information about idealized atmospheric variables at a given elevation.
   *)
  astro_atmosphere_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** Atmospheric pressure in pascals *)
    pressure: Double;
    (** Atmospheric temperature in kelvins *)
    temperature: Double;
    (** Atmospheric density relative to sea level *)
    density: Double;
  end;

  (**
   * @brief The result of a search for an astronomical event.
   *)
  astro_search_result_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The time at which a searched-for event occurs. *)
    time: astro_time_t;
  end;

  (**
   * @brief
   *      The dates and times of changes of season for a given calendar year.
   *      Call #Astronomy_Seasons to calculate this data structure for a given year.
   *)
  astro_seasons_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The date and time of the March equinox for the specified year. *)
    mar_equinox: astro_time_t;
    (** The date and time of the June soltice for the specified year. *)
    jun_solstice: astro_time_t;
    (** The date and time of the September equinox for the specified year. *)
    sep_equinox: astro_time_t;
    (** The date and time of the December solstice for the specified year. *)
    dec_solstice: astro_time_t;
  end;

  (**
   * @brief A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.
   *)
  astro_moon_quarter_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** 0=new moon, 1=first quarter, 2=full moon, 3=third quarter. *)
    quarter: Integer;
    (** The date and time of the lunar quarter. *)
    time: astro_time_t;
  end;

  (**
   * @brief A real value returned by a function whose ascending root is to be found.
   *
   * When calling #Astronomy_Search, the caller must pass in a callback function
   * compatible with the function-pointer type #astro_search_func_t
   * whose ascending root is to be found. That callback function must return astro_func_result_t.
   * If the function call is successful, it will set `status` to `ASTRO_SUCCESS` and `value`
   * to the numeric value appropriate for the given date and time.
   * If the call fails for some reason, it should set `status` to an appropriate error value
   * other than `ASTRO_SUCCESS`; in the error case, to guard against any possible misuse of `value`,
   * it is recommended to set `value` to `NAN`, though this is not strictly necessary.
   *)
  astro_func_result_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The value returned by a function whose ascending root is to be found. *)
    value: Double;
  end;

  (**
   * @brief A pointer to a function that is to be passed as a callback to #Astronomy_Search.
   *
   * The function #Astronomy_Search numerically solves for the time that a given event occurs.
   * An event is defined as the time when an arbitrary function transitions between having
   * a negative value and a non-negative value. This transition is called an *ascending root*.
   *
   * The type astro_search_func_t represents such a callback function that accepts a
   * custom `context` pointer and an astro_time_t representing the time to probe.
   * The function returns an astro_func_result_t that contains either a real
   * number in `value` or an error code in `status` that aborts the search.
   *
   * The `context` points to some data whose type varies depending on the callback function.
   * It can contain any auxiliary parameters (other than time) needed to evaluate the function.
   * For example, a function may pertain to a specific celestial body, in which case `context`
   * may point to a value of type astro_body_t. The `context` parameter is supplied by
   * the caller of #Astronomy_Search, which passes it along to every call to the callback function.
   * If the caller of `Astronomy_Search` knows that the callback function does not need a context,
   * it is safe to pass `NULL` as the context pointer.
   *)
  astro_search_func_t = function(context: Pointer; time: astro_time_t): astro_func_result_t; cdecl;

  (**
   * @brief A pointer to a function that calculates Delta T.
   *
   * Delta T is the discrepancy between times measured using an atomic clock
   * and times based on observations of the Earth's rotation, which is gradually
   * slowing down over time. Delta T = TT - UT, where
   * TT = Terrestrial Time, based on atomic time, and
   * UT = Universal Time, civil time based on the Earth's rotation.
   * Astronomy Engine defaults to using a Delta T function defined by
   * Espenak and Meeus in their "Five Millennium Canon of Solar Eclipses".
   * See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
   *)
  astro_deltat_func = function(ut: Double): Double; cdecl;

  (**
   * @brief Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.
   *)
  astro_visibility_t = (
    (** The body is best visible in the morning, before sunrise. *)
    VISIBLE_MORNING = 0,
    (** The body is best visible in the evening, after sunset. *)
    VISIBLE_EVENING = 1);
  Pastro_visibility_t = ^astro_visibility_t;

  (**
   * @brief
   *      Contains information about the visibility of a celestial body at a given date and time.
   *      See #Astronomy_Elongation for more detailed information about the members of this structure.
   *      See also #Astronomy_SearchMaxElongation for how to search for maximum elongation events.
   *)
  astro_elongation_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The date and time of the observation. *)
    time: astro_time_t;
    (** Whether the body is best seen in the morning or the evening. *)
    visibility: astro_visibility_t;
    (** The angle in degrees between the body and the Sun, as seen from the Earth. *)
    elongation: Double;
    (** The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth. *)
    ecliptic_separation: Double;
  end;

  (**
   * @brief Information about a celestial body crossing a specific hour angle.
   *
   * Returned by the function #Astronomy_SearchHourAngleEx to report information about
   * a celestial body crossing a certain hour angle as seen by a specified topocentric observer.
   *)
  astro_hour_angle_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The date and time when the body crosses the specified hour angle. *)
    time: astro_time_t;
    (** Apparent coordinates of the body at the time it crosses the specified hour angle. *)
    hor: astro_horizon_t;
  end;

  (**
   * @brief Information about the brightness and illuminated shape of a celestial body.
   *
   * Returned by the functions #Astronomy_Illumination and #Astronomy_SearchPeakMagnitude
   * to report the visual magnitude and illuminated fraction of a celestial body at a given date and time.
   *)
  astro_illum_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The date and time of the observation. *)
    time: astro_time_t;
    (** The visual magnitude of the body. Smaller values are brighter. *)
    mag: Double;
    (** The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth. *)
    phase_angle: Double;
    (** A value in the range [0.0, 1.0] indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth. *)
    phase_fraction: Double;
    (** The distance between the Sun and the body at the observation time. *)
    helio_dist: Double;
    (** For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0. *)
    ring_tilt: Double;
  end;

  (**
   * @brief The type of apsis: pericenter (closest approach) or apocenter (farthest distance).
   *)
  astro_apsis_kind_t = (
    (** The body is at its closest approach to the object it orbits. *)
    APSIS_PERICENTER = 0,
    (** The body is at its farthest distance from the object it orbits. *)
    APSIS_APOCENTER = 1,
    (** Undefined or invalid apsis. *)
    APSIS_INVALID = 2);
  Pastro_apsis_kind_t = ^astro_apsis_kind_t;

  (**
   * @brief An apsis event: pericenter (closest approach) or apocenter (farthest distance).
   *
   * For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
   * event where the orbiting body reaches its closest or farthest point from the primary body.
   * The closest approach is called *pericenter* and the farthest point is *apocenter*.
   *
   * More specific terminology is common for particular orbiting bodies.
   * The Moon's closest approach to the Earth is called *perigee* and its farthest
   * point is called *apogee*. The closest approach of a planet to the Sun is called
   * *perihelion* and the furthest point is called *aphelion*.
   *
   * This data structure is returned by #Astronomy_SearchLunarApsis and #Astronomy_NextLunarApsis
   * to iterate through consecutive alternating perigees and apogees.
   *)
  astro_apsis_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The date and time of the apsis. *)
    time: astro_time_t;
    (** Whether this is a pericenter or apocenter event. *)
    kind: astro_apsis_kind_t;
    (** The distance between the centers of the bodies in astronomical units. *)
    dist_au: Double;
    (** The distance between the centers of the bodies in kilometers. *)
    dist_km: Double;
  end;

  (**
   * @brief The different kinds of lunar/solar eclipses.
   *)
  astro_eclipse_kind_t = (
    (** No eclipse found. *)
    ECLIPSE_NONE = 0,
    (** A penumbral lunar eclipse. (Never used for a solar eclipse.) *)
    ECLIPSE_PENUMBRAL = 1,
    (** A partial lunar/solar eclipse. *)
    ECLIPSE_PARTIAL = 2,
    (** An annular solar eclipse. (Never used for a lunar eclipse.) *)
    ECLIPSE_ANNULAR = 3,
    (** A total lunar/solar eclipse. *)
    ECLIPSE_TOTAL = 4);
  Pastro_eclipse_kind_t = ^astro_eclipse_kind_t;

  (**
   * @brief Information about a lunar eclipse.
   *
   * Returned by #Astronomy_SearchLunarEclipse or #Astronomy_NextLunarEclipse
   * to report information about a lunar eclipse event.
   * If a lunar eclipse is found, `status` holds `ASTRO_SUCCESS` and the other fields are set.
   * If `status` holds any other value, it is an error code and the other fields are undefined.
   *
   * When a lunar eclipse is found, it is classified as penumbral, partial, or total.
   * Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed
   * by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
   * Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
   * Total eclipses occur when the entire Moon passes into the Earth's umbra.
   *
   * The `kind` field thus holds `ECLIPSE_PENUMBRAL`, `ECLIPSE_PARTIAL`, or `ECLIPSE_TOTAL`,
   * depending on the kind of lunar eclipse found.
   *
   * The `obscuration` field holds a value in the range [0, 1] that indicates what fraction
   * of the Moon's apparent disc area is covered by the Earth's umbra at the eclipse's peak.
   * This indicates how dark the peak eclipse appears. For penumbral eclipses, the obscuration
   * is 0, because the Moon does not pass through the Earth's umbra. For partial eclipses,
   * the obscuration is somewhere between 0 and 1. For total lunar eclipses, the obscuration is 1.
   *
   * Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
   *
   * Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
   * of the eclipse, which is half of the amount of time the eclipse spends in each
   * phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
   * By converting from minutes to days, and subtracting/adding with `center`, the caller
   * may determine the date and time of the beginning/end of each eclipse phase.
   *)
  astro_lunar_eclipse_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The type of lunar eclipse found. *)
    kind: astro_eclipse_kind_t;
    (** The peak fraction of the Moon's apparent disc that is covered by the Earth's umbra. *)
    obscuration: Double;
    (** The time of the eclipse at its peak. *)
    peak: astro_time_t;
    (** The semi-duration of the penumbral phase in minutes. *)
    sd_penum: Double;
    (** The semi-duration of the partial phase in minutes, or 0.0 if none. *)
    sd_partial: Double;
    (** The semi-duration of the total phase in minutes, or 0.0 if none. *)
    sd_total: Double;
  end;

  (**
   * @brief Reports the time and geographic location of the peak of a solar eclipse.
   *
   * Returned by #Astronomy_SearchGlobalSolarEclipse or #Astronomy_NextGlobalSolarEclipse
   * to report information about a solar eclipse event.
   * If a solar eclipse is found, `status` holds `ASTRO_SUCCESS` and `kind`, `peak`, and `distance`
   * have valid values. The `latitude` and `longitude` are set only for total and annular eclipses
   * (see more below).
   * If `status` holds any value other than `ASTRO_SUCCESS`, it is an error code;
   * in that case, `kind` holds `ECLIPSE_NONE` and all the other fields are undefined.
   *
   * The eclipse is classified as partial, annular, or total, depending on the
   * maximum amount of the Sun's disc obscured, as seen at the peak location
   * on the surface of the Earth.
   *
   * The `kind` field thus holds `ECLIPSE_PARTIAL`, `ECLIPSE_ANNULAR`, or `ECLIPSE_TOTAL`.
   * A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
   * An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
   * to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
   * A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
   * observes either a total or annular eclipse.
   *
   * If `kind` is `ECLIPSE_TOTAL` or `ECLIPSE_ANNULAR`, the `latitude` and `longitude`
   * fields give the geographic coordinates of the center of the Moon's shadow projected
   * onto the daytime side of the Earth at the instant of the eclipse's peak.
   * If `kind` has any other value, `latitude` and `longitude` are undefined and should
   * not be used.
   *
   * For total or annular eclipses, the `obscuration` field holds the fraction (0, 1]
   * of the Sun's apparent disc area that is blocked from view by the Moon's silhouette,
   * as seen by an observer located at the geographic coordinates `latitude`, `longitude`
   * at the darkest time `peak`. The value will always be 1 for total eclipses, and less than
   * 1 for annular eclipses.
   * For partial eclipses, `obscuration` is undefined and should not be used.
   * This is because there is little practical use for an obscuration value of
   * a partial eclipse without supplying a particular observation location.
   * Developers who wish to find an obscuration value for partial solar eclipses should therefore use
   * #Astronomy_SearchLocalSolarEclipse and provide the geographic coordinates of an observer.
   *)
  astro_global_solar_eclipse_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The type of solar eclipse found. *)
    kind: astro_eclipse_kind_t;
    (** The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only). *)
    obscuration: Double;
    (** The date and time when the solar eclipse is darkest. This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center. *)
    peak: astro_time_t;
    (** The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers. *)
    distance: Double;
    (** The geographic latitude at the center of the peak eclipse shadow. *)
    latitude: Double;
    (** The geographic longitude at the center of the peak eclipse shadow. *)
    longitude: Double;
  end;

  (**
   * @brief Holds a time and the observed altitude of the Sun at that time.
   *
   * When reporting a solar eclipse observed at a specific location on the Earth
   * (a "local" solar eclipse), a series of events occur. In addition
   * to the time of each event, it is important to know the altitude of the Sun,
   * because each event may be invisible to the observer if the Sun is below
   * the horizon.
   *
   * If `altitude` is negative, the event is theoretical only; it would be
   * visible if the Earth were transparent, but the observer cannot actually see it.
   * If `altitude` is positive but less than a few degrees, visibility will be impaired by
   * atmospheric interference (sunrise or sunset conditions).
   *)
  astro_eclipse_event_t = record
    (** The date and time of the event. *)
    time: astro_time_t;
    (** The angular altitude of the center of the Sun above/below the horizon, at `time`, corrected for atmospheric refraction and expressed in degrees. *)
    altitude: Double;
  end;

  (**
   * @brief Information about a solar eclipse as seen by an observer at a given time and geographic location.
   *
   * Returned by #Astronomy_SearchLocalSolarEclipse or #Astronomy_NextLocalSolarEclipse
   * to report information about a solar eclipse as seen at a given geographic location.
   * If a solar eclipse is found, `status` holds `ASTRO_SUCCESS` and the other fields are set.
   * If `status` holds any other value, it is an error code and the other fields are undefined.
   *
   * When a solar eclipse is found, it is classified as partial, annular, or total.
   * The `kind` field thus holds `ECLIPSE_PARTIAL`, `ECLIPSE_ANNULAR`, or `ECLIPSE_TOTAL`.
   * A partial solar eclipse is when the Moon does not line up directly enough with the Sun
   * to completely block the Sun's light from reaching the observer.
   * An annular eclipse occurs when the Moon's disc is completely visible against the Sun
   * but the Moon is too far away to completely block the Sun's light; this leaves the
   * Sun with a ring-like appearance.
   * A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
   * Sun just right to completely block all sunlight from reaching the observer.
   *
   * The `obscuration` field reports what fraction of the Sun's disc appears blocked
   * by the Moon when viewed by the observer at the peak eclipse time.
   * This is a value that ranges from 0 (no blockage) to 1 (total eclipse).
   * The obscuration value will be between 0 and 1 for partial eclipses and annular eclipses.
   * The value will be exactly 1 for total eclipses. Obscuration gives an indication
   * of how dark the eclipse appears.
   *
   * There are 5 "event" fields, each of which contains a time and a solar altitude.
   * Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
   * The fields `partial_begin` and `partial_end` are always set, and indicate when
   * the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
   * `total_begin` and `total_end` indicate when the total/annular phase begins/ends.
   * When an event field is valid, the caller must also check its `altitude` field to
   * see whether the Sun is above the horizon at that time. See #astro_eclipse_kind_t
   * for more information.
   *)
  astro_local_solar_eclipse_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The type of solar eclipse found: `ECLIPSE_PARTIAL`, `ECLIPSE_ANNULAR`, or `ECLIPSE_TOTAL`. *)
    kind: astro_eclipse_kind_t;
    (** The fraction of the Sun's apparent disc area obscured by the Moon at the eclipse peak. *)
    obscuration: Double;
    (** The time and Sun altitude at the beginning of the eclipse. *)
    partial_begin: astro_eclipse_event_t;
    (** If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise invalid. *)
    total_begin: astro_eclipse_event_t;
    (** The time and Sun altitude when the eclipse reaches its peak. *)
    peak: astro_eclipse_event_t;
    (** If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise invalid. *)
    total_end: astro_eclipse_event_t;
    (** The time and Sun altitude at the end of the eclipse. *)
    partial_end: astro_eclipse_event_t;
  end;

  (**
   * @brief Information about a transit of Mercury or Venus, as seen from the Earth.
   *
   * Returned by #Astronomy_SearchTransit or #Astronomy_NextTransit to report
   * information about a transit of Mercury or Venus.
   * A transit is when Mercury or Venus passes between the Sun and Earth so that
   * the other planet is seen in silhouette against the Sun.
   *
   * The `start` field reports the moment in time when the planet first becomes
   * visible against the Sun in its background.
   * The `peak` field reports when the planet is most aligned with the Sun,
   * as seen from the Earth.
   * The `finish` field reports the last moment when the planet is visible
   * against the Sun in its background.
   *
   * The calculations are performed from the point of view of a geocentric observer.
   *)
  astro_transit_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** Date and time at the beginning of the transit. *)
    start: astro_time_t;
    (** Date and time of the peak of the transit. *)
    peak: astro_time_t;
    (** Date and time at the end of the transit. *)
    finish: astro_time_t;
    (** Angular separation in arcminutes between the centers of the Sun and the planet at time `peak`. *)
    separation: Double;
  end;

  (**
   * @brief   Aberration calculation options.
   *
   * [Aberration](https://en.wikipedia.org/wiki/Aberration_of_light) is an effect
   * causing the apparent direction of an observed body to be shifted due to transverse
   * movement of the Earth with respect to the rays of light coming from that body.
   * This angular correction can be anywhere from 0 to about 20 arcseconds,
   * depending on the position of the observed body relative to the instantaneous
   * velocity vector of the Earth.
   *
   * Some Astronomy Engine functions allow optional correction for aberration by
   * passing in a value of this enumerated type.
   *
   * Aberration correction is useful to improve accuracy of coordinates of
   * apparent locations of bodies seen from the Earth.
   * However, because aberration affects not only the observed body (such as a planet)
   * but the surrounding stars, aberration may be unhelpful (for example)
   * for determining exactly when a planet crosses from one constellation to another.
   *)
  astro_aberration_t = (
    (** Request correction for aberration. *)
    ABERRATION = 0,
    (** Do not correct for aberration. *)
    NO_ABERRATION = 1);
  Pastro_aberration_t = ^astro_aberration_t;

  (**
   * @brief   Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.
   *
   * The Earth's equator is not always in the same plane due to precession and nutation.
   *
   * Sometimes it is useful to have a fixed plane of reference for equatorial coordinates
   * across different calendar dates.  In these cases, a fixed *epoch*, or reference time,
   * is helpful. Astronomy Engine provides the J2000 epoch for such cases.  This refers
   * to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.
   *
   * For some other purposes, it is more helpful to represent coordinates using the Earth's
   * equator exactly as it is on that date. For example, when calculating rise/set times
   * or horizontal coordinates, it is most accurate to use the orientation of the Earth's
   * equator at that same date and time. For these uses, Astronomy Engine allows *of-date*
   * calculations.
   *)
  astro_equator_date_t = (
    (** Represent equatorial coordinates in the J2000 epoch. *)
    EQUATOR_J2000 = 0,
    (** Represent equatorial coordinates using the Earth's equator at the given date and time. *)
    EQUATOR_OF_DATE = 1);
  Pastro_equator_date_t = ^astro_equator_date_t;

  (**
   * @brief Selects whether to search for a rise time or a set time.
   *
   * The #Astronomy_SearchRiseSetEx function finds the rise or set time of a body
   * depending on the value of its `direction` parameter.
   *)
  astro_direction_t = (
    (** Search for the time a body begins to rise above the horizon. *)
    DIRECTION_RISE = 1,
    (** Search for the time a body finishes sinking below the horizon. *)
    DIRECTION_SET = -1);
  Pastro_direction_t = ^astro_direction_t;

  (**
   * @brief Reports the constellation that a given celestial point lies within.
   *
   * The #Astronomy_Constellation function returns this struct
   * to report which constellation corresponds with a given point in the sky.
   * Constellations are defined with respect to the B1875 equatorial system
   * per IAU standard. Although `Astronomy.Constellation` requires J2000 equatorial
   * coordinates, the struct contains converted B1875 coordinates for reference.
   *)
  astro_constellation_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** 3-character mnemonic symbol for the constellation, e.g. "Ori". *)
    symbol: PUTF8Char;
    (** Full name of constellation, e.g. "Orion". *)
    name: PUTF8Char;
    (** Right ascension expressed in B1875 coordinates. *)
    ra_1875: Double;
    (** Declination expressed in B1875 coordinates. *)
    dec_1875: Double;
  end;

  (**
   * @brief Selects the output format of the function #Astronomy_FormatTime.
   *)
  astro_time_format_t = (
    (** Truncate to UTC calendar date only, e.g. `2020-12-31`. Buffer size must be at least 11 characters. *)
    TIME_FORMAT_DAY = 0,
    (** Round to nearest UTC minute, e.g. `2020-12-31T15:47Z`. Buffer size must be at least 18 characters. *)
    TIME_FORMAT_MINUTE = 1,
    (** Round to nearest UTC second, e.g. `2020-12-31T15:47:32Z`. Buffer size must be at least 21 characters. *)
    TIME_FORMAT_SECOND = 2,
    (** Round to nearest UTC millisecond, e.g. `2020-12-31T15:47:32.397Z`. Buffer size must be at least 25 characters. *)
    TIME_FORMAT_MILLI = 3);
  Pastro_time_format_t = ^astro_time_format_t;

  (**
   * @brief Lunar libration angles, returned by #Astronomy_Libration.
   *)
  astro_libration_t = record
    (** Sub-Earth libration ecliptic latitude angle, in degrees. *)
    elat: Double;
    (** Sub-Earth libration ecliptic longitude angle, in degrees. *)
    elon: Double;
    (** Moon's geocentric ecliptic latitude, in degrees. *)
    mlat: Double;
    (** Moon's geocentric ecliptic longitude, in degrees. *)
    mlon: Double;
    (** Distance between the centers of the Earth and Moon in kilometers. *)
    dist_km: Double;
    (** The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth. *)
    diam_deg: Double;
  end;

  (**
   * @brief Information about a body's rotation axis at a given time.
   *
   * This structure is returned by #Astronomy_RotationAxis to report
   * the orientation of a body's rotation axis at a given moment in time.
   * The axis is specified by the direction in space that the body's north pole
   * points, using angular equatorial coordinates in the J2000 system (EQJ).
   *
   * Thus `ra` is the right ascension, and `dec` is the declination, of the
   * body's north pole vector at the given moment in time. The north pole
   * of a body is defined as the pole that lies on the north side of the
   * [Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane),
   * regardless of the body's direction of rotation.
   *
   * The `spin` field indicates the angular position of a prime meridian
   * arbitrarily recommended for the body by the International Astronomical
   * Union (IAU).
   *
   * The fields `ra`, `dec`, and `spin` correspond to the variables
   * О±0, Оґ0, and W, respectively, from
   * [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
   *
   * The field `north` is a unit vector pointing in the direction of the body's north pole.
   * It is expressed in the equatorial J2000 system (EQJ).
   *)
  astro_axis_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The J2000 right ascension of the body's north pole direction, in sidereal hours. *)
    ra: Double;
    (** The J2000 declination of the body's north pole direction, in degrees. *)
    dec: Double;
    (** Rotation angle of the body's prime meridian, in degrees. *)
    spin: Double;
    (** A J2000 dimensionless unit vector pointing in the direction of the body's north pole. *)
    north: astro_vector_t;
  end;

  (**
   * @brief Holds the positions and velocities of Jupiter's major 4 moons.
   *
   * The #Astronomy_JupiterMoons function returns this struct
   * to report position and velocity vectors for Jupiter's largest 4 moons
   * Io, Europa, Ganymede, and Callisto. Each position vector is relative
   * to the center of Jupiter. Both position and velocity are oriented in
   * the EQJ system (that is, using Earth's equator at the J2000 epoch.)
   * The positions are expressed in astronomical units (AU),
   * and the velocities in AU/day.
   *)
  astro_jupiter_moons_t = record
    (** Jovicentric position and velocity of Io. *)
    io: astro_state_vector_t;
    (** Jovicentric position and velocity of Europa. *)
    europa: astro_state_vector_t;
    (** Jovicentric position and velocity of Ganymede. *)
    ganymede: astro_state_vector_t;
    (** Jovicentric position and velocity of Callisto. *)
    callisto: astro_state_vector_t;
  end;

  (**
   * @brief  Indicates whether a crossing through the ecliptic plane is ascending or descending.
   *)
  astro_node_kind_t = (
    (** Placeholder value for a missing or invalid node. *)
    INVALID_NODE = 0,
    (** The body passes through the ecliptic plane from south to north. *)
    ASCENDING_NODE = 1,
    (** The body passes through the ecliptic plane from north to south. *)
    DESCENDING_NODE = -1);
  Pastro_node_kind_t = ^astro_node_kind_t;

  (**
   * @brief Information about an ascending or descending node of a body.
   *
   * This structure is returned by #Astronomy_SearchMoonNode and #Astronomy_NextMoonNode
   * to report information about the center of the Moon passing through the ecliptic plane.
   *)
  astro_node_event_t = record
    (** `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. *)
    status: astro_status_t;
    (** The time when the body passes through the ecliptic plane. *)
    time: astro_time_t;
    (** Either `ASCENDING_NODE` or `DESCENDING_NODE`, depending on the direction of the ecliptic plane crossing. *)
    kind: astro_node_kind_t;
  end;

  Pastro_grav_sim_t = Pointer;
  PPastro_grav_sim_t = ^Pastro_grav_sim_t;

  (**
   * @brief A function for which to solve a light-travel time problem.
   *
   * The function #Astronomy_CorrectLightTravel solves a generalized
   * problem of deducing how far in the past light must have left
   * a target object to be seen by an observer at a specified time.
   * This function pointer type expresses an arbitrary position vector
   * as function of time. Such a function must be passed to
   * `Astronomy_CorrectLightTravel`.
   *)
  astro_position_func_t = function(context: Pointer; time: astro_time_t): astro_vector_t; cdecl;

function Astronomy_DeltaT_EspenakMeeus(ut: Double): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_DeltaT_EspenakMeeus';

function Astronomy_DeltaT_JplHorizons(ut: Double): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_DeltaT_JplHorizons';

procedure Astronomy_SetDeltaTFunction(func: astro_deltat_func); cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SetDeltaTFunction';

procedure Astronomy_Reset(); cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Reset';

function Astronomy_VectorLength(vector: astro_vector_t): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_VectorLength';

function Astronomy_AngleBetween(a: astro_vector_t; b: astro_vector_t): astro_angle_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_AngleBetween';

function Astronomy_BodyName(body: astro_body_t): PUTF8Char; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_BodyName';

function Astronomy_BodyCode(const name: PUTF8Char): astro_body_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_BodyCode';

function Astronomy_MakeObserver(latitude: Double; longitude: Double; height: Double): astro_observer_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_MakeObserver';

function Astronomy_CurrentTime(): astro_time_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_CurrentTime';

function Astronomy_MakeTime(year: Integer; month: Integer; day: Integer; hour: Integer; minute: Integer; second: Double): astro_time_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_MakeTime';

function Astronomy_TimeFromUtc(utc: astro_utc_t): astro_time_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_TimeFromUtc';

function Astronomy_UtcFromTime(time: astro_time_t): astro_utc_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_UtcFromTime';

function Astronomy_FormatTime(time: astro_time_t; format: astro_time_format_t; text: PUTF8Char; size: NativeUInt): astro_status_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_FormatTime';

function Astronomy_TimeFromDays(ut: Double): astro_time_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_TimeFromDays';

function Astronomy_TerrestrialTime(tt: Double): astro_time_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_TerrestrialTime';

function Astronomy_AddDays(time: astro_time_t; days: Double): astro_time_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_AddDays';

function Astronomy_SiderealTime(time: Pastro_time_t): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SiderealTime';

function Astronomy_HelioDistance(body: astro_body_t; time: astro_time_t): astro_func_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_HelioDistance';

function Astronomy_HelioVector(body: astro_body_t; time: astro_time_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_HelioVector';

function Astronomy_GeoVector(body: astro_body_t; time: astro_time_t; aberration: astro_aberration_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GeoVector';

function Astronomy_GeoMoon(time: astro_time_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GeoMoon';

function Astronomy_EclipticGeoMoon(time: astro_time_t): astro_spherical_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_EclipticGeoMoon';

function Astronomy_GeoMoonState(time: astro_time_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GeoMoonState';

function Astronomy_GeoEmbState(time: astro_time_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GeoEmbState';

function Astronomy_Libration(time: astro_time_t): astro_libration_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Libration';

function Astronomy_BaryState(body: astro_body_t; time: astro_time_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_BaryState';

function Astronomy_HelioState(body: astro_body_t; time: astro_time_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_HelioState';

function Astronomy_MassProduct(body: astro_body_t): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_MassProduct';

function Astronomy_PlanetOrbitalPeriod(body: astro_body_t): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_PlanetOrbitalPeriod';

function Astronomy_LagrangePoint(point: Integer; time: astro_time_t; major_body: astro_body_t; minor_body: astro_body_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_LagrangePoint';

function Astronomy_LagrangePointFast(point: Integer; major_state: astro_state_vector_t; major_mass: Double; minor_state: astro_state_vector_t; minor_mass: Double): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_LagrangePointFast';

function Astronomy_JupiterMoons(time: astro_time_t): astro_jupiter_moons_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_JupiterMoons';

function Astronomy_Equator(body: astro_body_t; time: Pastro_time_t; observer: astro_observer_t; equdate: astro_equator_date_t; aberration: astro_aberration_t): astro_equatorial_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Equator';

function Astronomy_ObserverVector(time: Pastro_time_t; observer: astro_observer_t; equdate: astro_equator_date_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_ObserverVector';

function Astronomy_ObserverState(time: Pastro_time_t; observer: astro_observer_t; equdate: astro_equator_date_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_ObserverState';

function Astronomy_VectorObserver(vector: Pastro_vector_t; equdate: astro_equator_date_t): astro_observer_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_VectorObserver';

function Astronomy_ObserverGravity(latitude: Double; height: Double): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_ObserverGravity';

function Astronomy_SunPosition(time: astro_time_t): astro_ecliptic_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SunPosition';

function Astronomy_Ecliptic(eqj: astro_vector_t): astro_ecliptic_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Ecliptic';

function Astronomy_EclipticLongitude(body: astro_body_t; time: astro_time_t): astro_angle_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_EclipticLongitude';

function Astronomy_Horizon(time: Pastro_time_t; observer: astro_observer_t; ra: Double; dec: Double; refraction: astro_refraction_t): astro_horizon_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Horizon';

function Astronomy_AngleFromSun(body: astro_body_t; time: astro_time_t): astro_angle_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_AngleFromSun';

function Astronomy_Elongation(body: astro_body_t; time: astro_time_t): astro_elongation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Elongation';

function Astronomy_SearchMaxElongation(body: astro_body_t; startTime: astro_time_t): astro_elongation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchMaxElongation';

function Astronomy_PairLongitude(body1: astro_body_t; body2: astro_body_t; time: astro_time_t): astro_angle_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_PairLongitude';

(** @endcond *)
function Astronomy_SearchRelativeLongitude(body: astro_body_t; targetRelLon: Double; startTime: astro_time_t): astro_search_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchRelativeLongitude';

function Astronomy_MoonPhase(time: astro_time_t): astro_angle_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_MoonPhase';

function Astronomy_SearchMoonPhase(targetLon: Double; startTime: astro_time_t; limitDays: Double): astro_search_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchMoonPhase';

function Astronomy_SearchMoonQuarter(startTime: astro_time_t): astro_moon_quarter_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchMoonQuarter';

function Astronomy_NextMoonQuarter(mq: astro_moon_quarter_t): astro_moon_quarter_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextMoonQuarter';

function Astronomy_SearchLunarEclipse(startTime: astro_time_t): astro_lunar_eclipse_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchLunarEclipse';

function Astronomy_NextLunarEclipse(prevEclipseTime: astro_time_t): astro_lunar_eclipse_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextLunarEclipse';

function Astronomy_SearchGlobalSolarEclipse(startTime: astro_time_t): astro_global_solar_eclipse_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchGlobalSolarEclipse';

function Astronomy_NextGlobalSolarEclipse(prevEclipseTime: astro_time_t): astro_global_solar_eclipse_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextGlobalSolarEclipse';

function Astronomy_SearchLocalSolarEclipse(startTime: astro_time_t; observer: astro_observer_t): astro_local_solar_eclipse_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchLocalSolarEclipse';

function Astronomy_NextLocalSolarEclipse(prevEclipseTime: astro_time_t; observer: astro_observer_t): astro_local_solar_eclipse_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextLocalSolarEclipse';

function Astronomy_SearchTransit(body: astro_body_t; startTime: astro_time_t): astro_transit_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchTransit';

function Astronomy_NextTransit(body: astro_body_t; prevTransitTime: astro_time_t): astro_transit_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextTransit';

function Astronomy_SearchMoonNode(startTime: astro_time_t): astro_node_event_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchMoonNode';

function Astronomy_NextMoonNode(prevNode: astro_node_event_t): astro_node_event_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextMoonNode';

function Astronomy_Search(func: astro_search_func_t; context: Pointer; t1: astro_time_t; t2: astro_time_t; dt_tolerance_seconds: Double): astro_search_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Search';

function Astronomy_SearchSunLongitude(targetLon: Double; startTime: astro_time_t; limitDays: Double): astro_search_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchSunLongitude';

function Astronomy_SearchHourAngleEx(body: astro_body_t; observer: astro_observer_t; hourAngle: Double; startTime: astro_time_t; direction: Integer): astro_hour_angle_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchHourAngleEx';

function Astronomy_HourAngle(body: astro_body_t; time: Pastro_time_t; observer: astro_observer_t): astro_func_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_HourAngle';

(** @endcond *)
function Astronomy_SearchRiseSetEx(body: astro_body_t; observer: astro_observer_t; direction: astro_direction_t; startTime: astro_time_t; limitDays: Double; metersAboveGround: Double): astro_search_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchRiseSetEx';

function Astronomy_SearchAltitude(body: astro_body_t; observer: astro_observer_t; direction: astro_direction_t; startTime: astro_time_t; limitDays: Double; altitude: Double): astro_search_result_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchAltitude';

function Astronomy_Atmosphere(elevationMeters: Double): astro_atmosphere_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Atmosphere';

function Astronomy_RotationAxis(body: astro_body_t; time: Pastro_time_t): astro_axis_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_RotationAxis';

function Astronomy_Seasons(year: Integer): astro_seasons_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Seasons';

function Astronomy_Illumination(body: astro_body_t; time: astro_time_t): astro_illum_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Illumination';

function Astronomy_SearchPeakMagnitude(body: astro_body_t; startTime: astro_time_t): astro_illum_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchPeakMagnitude';

function Astronomy_SearchLunarApsis(startTime: astro_time_t): astro_apsis_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchLunarApsis';

function Astronomy_NextLunarApsis(apsis: astro_apsis_t): astro_apsis_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextLunarApsis';

function Astronomy_SearchPlanetApsis(body: astro_body_t; startTime: astro_time_t): astro_apsis_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SearchPlanetApsis';

function Astronomy_NextPlanetApsis(body: astro_body_t; apsis: astro_apsis_t): astro_apsis_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_NextPlanetApsis';

function Astronomy_IdentityMatrix(): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_IdentityMatrix';

function Astronomy_InverseRotation(rotation: astro_rotation_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_InverseRotation';

function Astronomy_CombineRotation(a: astro_rotation_t; b: astro_rotation_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_CombineRotation';

function Astronomy_Pivot(rotation: astro_rotation_t; axis: Integer; angle: Double): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Pivot';

function Astronomy_VectorFromSphere(sphere: astro_spherical_t; time: astro_time_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_VectorFromSphere';

function Astronomy_SphereFromVector(vector: astro_vector_t): astro_spherical_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_SphereFromVector';

function Astronomy_EquatorFromVector(vector: astro_vector_t): astro_equatorial_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_EquatorFromVector';

function Astronomy_VectorFromHorizon(sphere: astro_spherical_t; time: astro_time_t; refraction: astro_refraction_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_VectorFromHorizon';

function Astronomy_HorizonFromVector(vector: astro_vector_t; refraction: astro_refraction_t): astro_spherical_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_HorizonFromVector';

function Astronomy_RotateVector(rotation: astro_rotation_t; vector: astro_vector_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_RotateVector';

function Astronomy_RotateState(rotation: astro_rotation_t; state: astro_state_vector_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_RotateState';

function Astronomy_Rotation_EQD_EQJ(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQD_EQJ';

function Astronomy_Rotation_EQD_ECL(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQD_ECL';

function Astronomy_Rotation_EQD_ECT(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQD_ECT';

function Astronomy_Rotation_EQD_HOR(time: Pastro_time_t; observer: astro_observer_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQD_HOR';

function Astronomy_Rotation_EQJ_EQD(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQJ_EQD';

function Astronomy_Rotation_EQJ_ECT(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQJ_ECT';

function Astronomy_Rotation_EQJ_ECL(): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQJ_ECL';

function Astronomy_Rotation_EQJ_HOR(time: Pastro_time_t; observer: astro_observer_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQJ_HOR';

function Astronomy_Rotation_ECL_EQD(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_ECL_EQD';

function Astronomy_Rotation_ECL_EQJ(): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_ECL_EQJ';

function Astronomy_Rotation_ECL_HOR(time: Pastro_time_t; observer: astro_observer_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_ECL_HOR';

function Astronomy_Rotation_ECT_EQJ(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_ECT_EQJ';

function Astronomy_Rotation_ECT_EQD(time: Pastro_time_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_ECT_EQD';

function Astronomy_Rotation_HOR_EQD(time: Pastro_time_t; observer: astro_observer_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_HOR_EQD';

function Astronomy_Rotation_HOR_EQJ(time: Pastro_time_t; observer: astro_observer_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_HOR_EQJ';

function Astronomy_Rotation_HOR_ECL(time: Pastro_time_t; observer: astro_observer_t): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_HOR_ECL';

function Astronomy_Rotation_EQJ_GAL(): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_EQJ_GAL';

function Astronomy_Rotation_GAL_EQJ(): astro_rotation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Rotation_GAL_EQJ';

function Astronomy_Refraction(refraction: astro_refraction_t; altitude: Double): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Refraction';

function Astronomy_InverseRefraction(refraction: astro_refraction_t; bent_altitude: Double): Double; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_InverseRefraction';

function Astronomy_Constellation(ra: Double; dec: Double): astro_constellation_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_Constellation';

function Astronomy_GravSimInit(simOut: PPastro_grav_sim_t; originBody: astro_body_t; time: astro_time_t; numBodies: Integer; const bodyStateArray: Pastro_state_vector_t): astro_status_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimInit';

function Astronomy_GravSimUpdate(sim: Pastro_grav_sim_t; time: astro_time_t; numBodies: Integer; bodyStateArray: Pastro_state_vector_t): astro_status_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimUpdate';

function Astronomy_GravSimBodyState(sim: Pastro_grav_sim_t; body: astro_body_t): astro_state_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimBodyState';

function Astronomy_GravSimTime(sim: Pastro_grav_sim_t): astro_time_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimTime';

function Astronomy_GravSimNumBodies(sim: Pastro_grav_sim_t): Integer; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimNumBodies';

function Astronomy_GravSimOrigin(sim: Pastro_grav_sim_t): astro_body_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimOrigin';

procedure Astronomy_GravSimSwap(sim: Pastro_grav_sim_t); cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimSwap';

procedure Astronomy_GravSimFree(sim: Pastro_grav_sim_t); cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_GravSimFree';

(**
 * @brief A function for which to solve a light-travel time problem.
 *
 * The function #Astronomy_CorrectLightTravel solves a generalized
 * problem of deducing how far in the past light must have left
 * a target object to be seen by an observer at a specified time.
 * This function pointer type expresses an arbitrary position vector
 * as function of time. Such a function must be passed to
 * `Astronomy_CorrectLightTravel`.
 *)
 
function Astronomy_CorrectLightTravel(context: Pointer; func: astro_position_func_t; time: astro_time_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_CorrectLightTravel';

function Astronomy_BackdatePosition(time: astro_time_t; observerBody: astro_body_t; targetBody: astro_body_t; aberration: astro_aberration_t): astro_vector_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_BackdatePosition';

function Astronomy_DefineStar(body: astro_body_t; ra: Double; dec: Double; distanceLightYears: Double): astro_status_t; cdecl;
  external LIB_ASTRONOMY name _PU + 'Astronomy_DefineStar';

implementation

end. // astronomy unit
