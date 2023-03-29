(* --------------------------------------------------------------------------- *)
(* EclTimer *)
(* Computation of local circumstances for solar eclipses *)
(* References: *)
(* Explanatory Supplement to the Astronomical Ephemeris and the *)
(* American Ephemeris and Nautical Almanac; Her Majesty's Stationery *)
(* Office (1961). *)
(* P.K.Seidelmann; Explanatory Supplement to the  Astronomical Almananc; *)
(* University Science Press, Mill Valley, California (1992). *)
(* --------------------------------------------------------------------------- *)

program EclTimer(Input, Output);

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.Spheric,
  Apc.Sun,
  Apc.Moon,
  Apc.Time;

const
  MAX_TP_DEG = 8; (* Degree of Tschebyscheff polynomials *)
  TDEVELOP = 0.5E-4; (* Approx. 2d in julian centuries *)
  K_MOON = 0.2725076; (* Ratio moon/earth radii *)
  K_SUN = 109.1227; (* Ratio sun/earth radii *)
  AE = 23454.77992; (* 1AU in earth radii *)

type
  SHADOW_TYPE = (PENUMBRA, UMBRA);
  PHASE_TYPE = (NO_ECLIPSE, PARTIAL, ANNULAR, TOTAL);
  TPINDX = (RAM, DEM, RM, RAS, DES, RS);
  TPOLYEPH = array [TPINDX] of TPolynomCheb;
  REAL4 = array [1 .. 4] of Double;

  Tfunction = function(X: Double): Double; (* function prototype y=f(x) *)
  (* for passing a function as *)
  (* a subroutine parameter *)
  (* in Turbo Pascal *)

var
  TNEWMOON, ETDIFUT: Double;
  LAMBDA, PHI, RCOSPHI, RSINPHI: Double;
  CONTACT_TIMES, POSANG_P, POSANG_V: REAL4;
  T_MAX, MAX_MAGN, MAX_OBSC: Double;
  MAX_PHASE: PHASE_TYPE;
  POLYEPH: TPOLYEPH;

  (* --------------------------------------------------------------------------- *)
  (* *)
  (* MOD360: reduces the (angular) argument X to the interval [0..360] *)
  (* *)
  (* --------------------------------------------------------------------------- *)

function MOD360(X: Double): Double;
begin
  if (X > 0.0) then
    X := X - 360.0 * Trunc(X / 360.0)
  else
    X := X - 360.0 * Trunc(X / 360.0) + 360.0;
  MOD360 := X;
end;

(* --------------------------------------------------------------------------- *)
(* ETminUT: *)
(* Difference between ephemeris time and universal time (polynomial *)
(* representation for the years 1900-1995 and approximation for times *)
(* before 1900. *)
(* T      Time in julian centuries since J2000 (=(JD-2451545)/36525) *)
(* DTSEC  DT=ET-UT in sec (only if VALID=True) *)
(* VALID  True if T is before 1995 *)
(* Note: *)
(* Modified version of routine ETminUT from "Astronomy on the PC". *)
(* For historic periods the approximate relation given in *)
(* F.R.Stephenson, "Pre-Telescopic Astronomical Observations", *)
(* in "Tidal Friction and the Earth's Rotation", *)
(* ed. P.Brosche & J.Suendermann, Springer-Verlag 1978 *)
(* is used. *)
(* --------------------------------------------------------------------------- *)

procedure ETminUT(T: Double; var DTSEC: Double; var VALID: Boolean);

begin

  VALID := (T <= -0.05);
  if (VALID) then
  begin
    if (-1.0 <= T) then
      DTSEC := ((((-449.50 * T - 783.42) * T - 387.70) * T) + 13.34) * T + 62.14
    else
      DTSEC := (38.3 * T + 190.6) * T + 172.3;
  end;

end;

(* --------------------------------------------------------------------------- *)
(* Pegasus: Root finder using the Pegasus method *)
(* FUNC      : function to be examined of the form y=f(x) *)
(* LOWERBOUND: lower boundary for search interval *)
(* UPPERBOUND: upper boundary for search interval *)
(* ACCURACY  : desired accuracy for root bracketing *)
(* SUCCESS   : True, if a root has been found in the search interval *)
(* within a given maximum number of iterations *)
(* ROOT      : root of f(x), if SUCCESS=True *)
(* For the application of the method it is required that f(x) has different *)
(* sign at the interval boundaries. Otherwise the routine returns with *)
(* SUCCESS=False. *)
(* References: *)
(* Dowell M., Jarratt P., 'A modified Regula Falsi Method for Computing *)
(* the root of an equation', BIT 11, p.168-174 (1971). *)
(* Dowell M., Jarratt P., 'The "Pegasus Method for Computing the root *)
(* of an equation', BIT 12, p.503-508 (1972). *)
(* G.Engeln-Muellges, F.Reutter, 'Formelsammlung zur Numerischen *)
(* Mathematik mit FORTRAN77-Programmen', Bibliogr. Institut, *)
(* Zuerich (1986). *)
(* --------------------------------------------------------------------------- *)

procedure Pegasus(FUNC: Tfunction; LOWERBOUND: Double; UPPERBOUND: Double; ACCURACY: Double;
  var SUCCESS: Boolean; var ROOT: Double);

const
  MAXITERAT = 30;

var
  X1, X2, X3, F1, F2, F3: Double;
  COUNT: integer;

begin

  X1 := LOWERBOUND;
  F1 := FUNC(X1);
  X2 := UPPERBOUND;
  F2 := FUNC(X2);

  COUNT := 0;
  ROOT := X1;
  SUCCESS := False;

  if (F1 * F2 < 0.0) then
    repeat

      (* Approximation of the root by interpolation *)

      X3 := X2 - F2 / ((F2 - F1) / (X2 - X1));
      F3 := FUNC(X3);

      (* Replace (x1,f2) and (x2,f2) by new values, such that *)
      (* the root is again within the interval [x1,x2] *)

      if (F3 * F2 <= 0.0) then (* Root in [x2,x3] *)
      begin
        X1 := X2;
        F1 := F2; (* Replace (x1,f1) by (x2,f2) *)
        X2 := X3;
        F2 := F3; (* Replace (x2,f2) by (x3,f3) *)
      end
      else (* Root in [x1,x3] *)
      begin
        F1 := F1 * F2 / (F2 + F3); (* Replace (x1,f1) by (x1,f1') *)
        X2 := X3;
        F2 := F3; (* Replace (x2,f2) by (x3,f3) *)
      end;

      if (abs(F1) < abs(F2)) then
        ROOT := X1
      else
        ROOT := X2;

      SUCCESS := (abs(X2 - X1) <= ACCURACY);

      COUNT := COUNT + 1;

    until (SUCCESS or (COUNT = MAXITERAT));

end; (* Pegasus *)

(* --------------------------------------------------------------------------- *)
(* *)
(* GetInput: Input of the time of new moon and the observer's co-ordinates *)
(* *)
(* TNEWMOON: time of new moon in julian centuries since J2000 (UT) *)
(* ETDIFUT : difference between ephemeris time and universal time [s] *)
(* LAMBDA  : geographic longitude of the observer [deg] *)
(* PHI     : geographic latitude of the observer [deg] *)
(* *)
(* --------------------------------------------------------------------------- *)

procedure GetInput(var TNEWMOON, ETDIFUT, LAMBDA, PHI: Double);

var
  DAY, MONTH, YEAR: integer;
  UT: Double;
  DTVALID: Boolean;

begin

  writeln;
  writeln('     EclTimer: Local circumstances of solar eclipses');
  writeln('                     version 93/07/01               ');
  writeln('       (c) 1993 Thomas Pfleger, Oliver Montenbruck  ');
  writeln;

  write(' Date of new moon (YYYY MM DD UT)            ... ');
  Readln(YEAR, MONTH, DAY, UT);
  TNEWMOON := (MJD(DAY, MONTH, YEAR, UT) - 51544.5) / 36525.0;

  ETminUT(TNEWMOON, ETDIFUT, DTVALID);
  if (DTVALID) then
    write(' Difference ET-UT (proposal:', Trunc(ETDIFUT + 0.5):3, ' sec)         ... ')
  else
    write(' Difference ET-UT (sec)                      ... ');
  Readln(ETDIFUT);

  write(' Observer''s coordinates: longitude (>0 west) ... ');
  Readln(LAMBDA);
  write('                         latitude            ... ');
  Readln(PHI);

end; (* GetInput *)

(* --------------------------------------------------------------------------- *)
(* Bessel *)
(* computes the orientation of the fundamental plane, the coordinates *)
(* of the shadow center and the parameters of the shadow cone *)
(* T_UT     Time in julan centuries UT since J2000 *)
(* ETDIFUT  Difference between ephemeris time and universal time [s] *)
(* POLYEPH  Chebyshev polynomials of solar and lunar coordinates *)
(* IJK      Unit vectors of the fundamental plane (equatorial) *)
(* XSH, YSH Coordinates of the shadow center on the fundamental plane *)
(* (in earth radii) *)
(* F1       Half angle of the penumbra cone (in deg) *)
(* L1       Radius of the penumbra cone on the fundamental plane *)
(* (in earth radii) *)
(* F2       Half angle of the umbra cone (in deg) *)
(* L2       Radius of the umbra cone on the fundamental plane *)
(* (in earth radii) *)
(* --------------------------------------------------------------------------- *)

procedure Bessel(T_UT, ETDIFUT: Double; POLYEPH: TPOLYEPH; var IJK: Mat3X; var XSH, YSH: Double;
  var F1, L1, F2, L2: Double);

var
  T_ET, DIST, ZSH: Double;
  RA_SUN, DEC_SUN, R_SUN: Double;
  RA_MOON, DEC_MOON, R_MOON: Double;
  SUN_POS, MOON_POS, MS: Vector;
  I: Index;

begin

  (* Calculate cartesian coordinates of sun and moon *)
  T_ET := T_UT + ETDIFUT / (86400.0 * 36525.0);

  RA_SUN := T_Eval(POLYEPH[RAS], T_ET);
  DEC_SUN := T_Eval(POLYEPH[DES], T_ET);
  R_SUN := T_Eval(POLYEPH[RS], T_ET) * AE;
  Cart(R_SUN, DEC_SUN, RA_SUN, SUN_POS[X], SUN_POS[Y], SUN_POS[Z]);

  RA_MOON := T_Eval(POLYEPH[RAM], T_ET);
  DEC_MOON := T_Eval(POLYEPH[DEM], T_ET);
  R_MOON := T_Eval(POLYEPH[RM], T_ET);
  Cart(R_MOON, DEC_MOON, RA_MOON, MOON_POS[X], MOON_POS[Y], MOON_POS[Z]);

  (* The sun-moon direction unit vector becomes vector *)
  (* IJK[3,*] of the triade defining the fundamental plane *)
  for I := X to Z do
    MS[I] := SUN_POS[I] - MOON_POS[I];
  DIST := Norm(MS); (* Distance Sun-Moon *)
  for I := X to Z do
    IJK[3, I] := MS[I] / DIST;

  (* The unit vector I=IJK[1,*] lies in the equatorial *)
  (* plane and is perpendicular to K *)

  DIST := SQRT(IJK[3, X] * IJK[3, X] + IJK[3, Y] * IJK[3, Y]);

  IJK[1, X] := -IJK[3, Y] / DIST;
  IJK[1, Y] := +IJK[3, X] / DIST;
  IJK[1, Z] := 0.0;

  (* The unit vector J=IJK[2,*] is perpendicular to K and I *)
  Cross(IJK[3], IJK[1], IJK[2]);

  (* Lunar coordinates in the fundamental plane system *)
  XSH := Dot(MOON_POS, IJK[1]);
  YSH := Dot(MOON_POS, IJK[2]);
  ZSH := Dot(MOON_POS, IJK[3]);

  (* Sun-moon distance *)
  for I := X to Z do
    MS[I] := SUN_POS[I] - MOON_POS[I];
  DIST := Norm(MS);

  (* Shadow cones *)
  F1 := ASN((K_SUN + K_MOON) / DIST);
  L1 := ZSH * (K_SUN + K_MOON) / DIST + K_MOON;

  F2 := ASN((K_SUN - K_MOON) / DIST);
  L2 := ZSH * (K_SUN - K_MOON) / DIST - K_MOON;

end; (* Bessel *)

(* --------------------------------------------------------------------------- *)
(* Observer *)
(* projects the location of the observer onto the fundamental plane *)
(* T_UT                    Time (julian centuries UT since J2000) *)
(* LAMBDA,RCOSPHI,RSINPHI  Observer's geozentric coordinates *)
(* IJK                     Unit vectors of the fundamental plane at T_UT *)
(* XI,ETA,ZETA             observer's coordinates on the fundamental plane *)
(* --------------------------------------------------------------------------- *)

procedure Observer(T_UT: Double; LAMBDA, RCOSPHI, RSINPHI: Double; IJK: Mat3X;
  var XI, ETA, ZETA: Double);

var
  LHA_OBS: Double;
  ROBS: Vector;

begin

  (* Compute local hour angle of observer *)

  LHA_OBS := 15.0 * LMST(T_UT * 36525.0 + 51544.5, LAMBDA);

  (* Determine equatorial cartesian coordinates of the observer *)

  ROBS[X] := RCOSPHI * CS(LHA_OBS);
  ROBS[Y] := RCOSPHI * SN(LHA_OBS);
  ROBS[Z] := RSINPHI;

  (* Projection into the coordinate system of the fundamental plane *)

  XI := Dot(ROBS, IJK[1]);
  ETA := Dot(ROBS, IJK[2]);
  ZETA := Dot(ROBS, IJK[3]);

end; (* Observer *)

(* --------------------------------------------------------------------------- *)
(* *)
(* SHADOW_DIST: *)
(* *)
(* Computes the shadow distance function f(t)=D(t)**2-L(t)**2, where D is *)
(* the distance of the observer from the shadow axis, while L is the *)
(* radius of the shadow cone. For f(t)=0 the observer touches the shadow *)
(* cone (unit Earth radii squared). *)
(* *)
(* Shadow                 Umbra or penumbra *)
(* T_UT                   Time (julian centuries UT since J2000) *)
(* ETDIFUT                Difference dT=ET-UT (in seconds) *)
(* POLYEPH                Chebyshev polynomials of solar and lunar coord. *)
(* LAMBDA,RCOSPHI,RSINPHI Geocentric coordinates of the observer *)
(* *)
(* --------------------------------------------------------------------------- *)

function SHADOW_DIST(Shadow: SHADOW_TYPE; T_UT, ETDIFUT: Double; POLYEPH: TPOLYEPH;
  LAMBDA, RCOSPHI, RSINPHI: Double): Double;

var
  XI, ETA, ZETA: Double;
  XSH, YSH: Double;
  F1, L1, F2, L2, LL: Double;
  IJK: Mat3X;

begin

  (* Fundamental plane coordinates of the shadow and the observer *)

  Bessel(T_UT, ETDIFUT, POLYEPH, IJK, XSH, YSH, F1, L1, F2, L2);

  Observer(T_UT, LAMBDA, RCOSPHI, RSINPHI, IJK, XI, ETA, ZETA);

  (* Shadow radius at the observer's place (LL<0 for a total eclipse) *)

  case (Shadow) of
    PENUMBRA:
      LL := L1 - ZETA * TN(F1);
    UMBRA:
      LL := L2 - ZETA * TN(F2);
  end;

  (* Shadow distance function *)

  SHADOW_DIST := (XSH - XI) * (XSH - XI) + (YSH - ETA) * (YSH - ETA) - LL * LL;

end; (* SHADOW_DIST *)

(* --------------------------------------------------------------------------- *)
(* *)
(* UMBRAS_DIST: *)
(* *)
(* Computes the observer's distance from the umbra as a function of time *)
(* T_UT. Interface routine for use with Pegasus. *)
(* Global variables: ETDIFUT,POLYEPH,LAMBDA,RCOSPHI,RSINPHI. *)
(* --------------------------------------------------------------------------- *)
(*$F+ Turbo Pascal compilation directive: force far call *)

function UMBRAS_DIST(T_UT: Double): Double;
begin
  UMBRAS_DIST := SHADOW_DIST(UMBRA, T_UT, ETDIFUT, POLYEPH, LAMBDA, RCOSPHI, RSINPHI);
end;
(*$F- Turbo Pascal compilation directive: end far call *)
(* --------------------------------------------------------------------------- *)
(* *)
(* CONTACTS *)
(* *)
(* Determines the phase, magnitude and times of contacts of a solar *)
(* eclipse for a specified location of the observer based on the time *)
(* of new moon. *)
(* *)
(* T_NEWMOON                Approximate time of new moon *)
(* ETDIFUT                  Difference ephemeris time - universal time (s) *)
(* LAMBDA,RCOSPHI,RSINPHI   Geocentric coordinates of the observer *)
(* POLYEPH                  Chebychev polynomials for solar and lunar *)
(* coordinates *)
(* CONTACT_TIMES            Times of 1st to 4th contact (UT) *)
(* T_MAX                    Time of maximum eclipse (UT) *)
(* MAG                      Magnitude of the eclipse (eclipsed fraction *)
(* of the solar diameter at the maximum) *)
(* OBSCUR                   Degree of obscuration (area of sun = 1) *)
(* PHASE                    Phase at maximum eclipse *)
(* *)
(* All times in julian centuries since J2000. *)
(* *)
(* *)
(* Description: *)
(* *)
(* For computing the times of contact the roots of f(t)=D(t)**2-L(t)**2 are *)
(* determined, where D is the observer's distance from the shadow axis, *)
(* L is the radius of the shadow cone. For f(t)=0 the observer lies on the *)
(* surface of the shadow cone and sees the touching rims of sun and moon. *)
(* *)
(* Using quadratic interpolation, the times of 1st and 4th contact are first *)
(* computed, at which the observer touches the penumbra cone. Simultaneously *)
(* the time of maximum eclipse is obtained, at which the distance between *)
(* observer and shadow axis attains a minimum. *)
(* *)
(* If the eclipse is at least a partial one, CONTACT computes the magnitude *)
(* of the eclipse and the obscuration of the sun by the moon for the time of *)
(* maximum eclipse. *)
(* *)
(* In case the observer's distance from the shadow axis and the umbra cone *)
(* diameter at the time of maximum eclipse result in a total or annular *)
(* eclipse, the times of 2nd and 3rd contact may subsequently be derived. *)
(* Since the duration of the total or annular phase never exceeds a value of *)
(* approximately 13 minutes, the 2nd and 3rd contact take place in a *)
(* corresponding interval before or after the maximum. Making use of this *)
(* fact, one may use a method like the regula falsi or the Pegasus method, *)
(* which converge more rapidly than quadratic interpolation, for computing *)
(* the contact times. *)
(* *)
(* --------------------------------------------------------------------------- *)

procedure CONTACTS(T_NEWMOON, ETDIFUT: Double; LAMBDA, RCOSPHI, RSINPHI: Double; POLYEPH: TPOLYEPH;
  var CONTACT_TIMES: REAL4; var T_MAX, MAG, OBSCUR: Double; var PHASE: PHASE_TYPE);

const
  RANGE = 5.0; (* Search interval = +/-(RANGE+DTAB) [h] *)
  DTAB = 0.10; (* Step size [h] *)
  CENT = 876600.0; (* 24*36525 hours per century *)
  EPS = 1.0E-10; (* Accuracy for contact times (approx. 0.3s) *)


var
  I, N_ROOTS, N_CONTACTS: integer;
  MOON_POS: Vector;
  IJK: Mat3X;
  T_UT, DT: Double;
  PD_MINUS, PD_0, PD_PLUS, XE, YE: Double;
  ROOT1, ROOT2: Double;
  X1, Y1, XI, ETA, ZETA: Double;
  F1, F2, L1, L2, LL1, LL2, M: Double;
  RA_MOON, DEC_MOON, R_MOON: Double;
  A, B, C, S: Double;
  OK: Boolean;

begin (* CONTACTS *)

  (* Initialisation *)

  for I := 1 to 4 do
    CONTACT_TIMES[I] := T_NEWMOON;

  (* ----------------------------------------------------------------------- *)
  (* *)
  (* Search for 1st and 4th contact by quadratic interpolation. *)
  (* *)
  (* The variation of the observer's distance from the penumbra in the *)
  (* fundamental plane is represented by three equidistant function values *)
  (* PD_MINUS, PD_0 and PD_PLUS separated by a timespan of DTAB. *)
  (* An interpolating parabola through these points and its roots is then *)
  (* computed, which yields the times of the 1st and 4th contact with *)
  (* moderate effort, since the interpolated function (i.e. the square of *)
  (* the shadow distance) is sufficiently smooth. *)
  (* *)
  (* ----------------------------------------------------------------------- *)

  N_CONTACTS := 0;

  T_UT := T_NEWMOON + (-RANGE - DTAB) / CENT;
  PD_PLUS := SHADOW_DIST(PENUMBRA, T_UT, ETDIFUT, POLYEPH, LAMBDA, RCOSPHI, RSINPHI);

  T_UT := T_NEWMOON + (-RANGE - 2.0 * DTAB) / CENT;

  repeat

    (* Compute next time step *)

    T_UT := T_UT + 2.0 * DTAB / CENT;

    (* Compute square of the shadow distance at times T_UT-DTAB, *)
    (* T_UT and T_UT+DTAB and interpolate *)

    PD_MINUS := PD_PLUS;
    PD_0 := SHADOW_DIST(PENUMBRA, T_UT, ETDIFUT, POLYEPH, LAMBDA, RCOSPHI, RSINPHI);
    PD_PLUS := SHADOW_DIST(PENUMBRA, T_UT + DTAB / CENT, ETDIFUT, POLYEPH, LAMBDA, RCOSPHI,
      RSINPHI);

    Quad(PD_MINUS, PD_0, PD_PLUS, XE, YE, ROOT1, ROOT2, N_ROOTS);

    (* Store number of contacts found so far and compute contact times *)

    N_CONTACTS := N_CONTACTS + N_ROOTS;

    case (N_ROOTS) of
      1:
        if (N_CONTACTS = 1) then
          CONTACT_TIMES[1] := T_UT + (DTAB * ROOT1) / CENT (* 1st contact *)
        else
          CONTACT_TIMES[4] := T_UT + (DTAB * ROOT1) / CENT; (* 4th contact *)
      2:
        begin
          CONTACT_TIMES[1] := T_UT + (DTAB * ROOT1) / CENT; (* 1st contact *)
          CONTACT_TIMES[4] := T_UT + (DTAB * ROOT2) / CENT; (* 4th contact *)
        end;
    end;

    (* If the minimum of the shadow distance lies within the present *)
    (* interval, a sufficiently accurate value of the time of maximum *)
    (* eclipse is determined from its position without further iteration. *)

    if (-1.0 < XE) and (XE < 1.0) then
    begin
      T_MAX := T_UT + (DTAB * XE) / CENT;
    end;

  until ((T_NEWMOON + RANGE / CENT <= T_UT) or (N_CONTACTS = 2));

  (* Determine type of eclipse *)

  if (N_CONTACTS = 0) then
    PHASE := NO_ECLIPSE (* No eclipse found *)
  else
    PHASE := PARTIAL; (* Eclipse is at least partial *)

  (* ----------------------------------------------------------------------- *)
  (* *)
  (* Degree of obscuration and magnitude of the eclipse *)
  (* *)
  (* ----------------------------------------------------------------------- *)

  if (PHASE > NO_ECLIPSE) then
  begin

    (* Coordinates of shadow and observer on the fundamental plane *)
    (* at maximum eclipse *)

    Bessel(T_MAX, ETDIFUT, POLYEPH, IJK, X1, Y1, F1, L1, F2, L2);
    Observer(T_MAX, LAMBDA, RCOSPHI, RSINPHI, IJK, XI, ETA, ZETA);

    (* Distance between observer and shadow axis *)

    M := SQRT((XI - X1) * (XI - X1) + (ETA - Y1) * (ETA - Y1));

    (* Penumbra and umbra radius at the observer's place *)
    (* (L2<0 for a total eclipse!) *)

    LL1 := L1 - ZETA * TN(F1);
    LL2 := L2 - ZETA * TN(F2);

    (* Eclipse type *)

    if (M < +LL2) then
      PHASE := ANNULAR;
    if (M < -LL2) then
      PHASE := TOTAL;

    (* Eclipse magnitude *)

    if (PHASE = PARTIAL) then
      MAG := (LL1 - M) / (LL1 + LL2) (* Penumbra region *)
    else
      MAG := (LL1 - LL2) / (LL1 + LL2); (* Umbra region *)

    (* Compute degree of obscuration *)

    case (PHASE) of
      NO_ECLIPSE:
        OBSCUR := 0.0;
      PARTIAL:
        begin
          B := (Pi / 180.0) * ACS((LL1 * LL2 + M * M) / (M * (LL1 + LL2)));
          C := (Pi / 180.0) * ACS((LL1 * LL1 + LL2 * LL2 - 2 * M * M) / (LL1 * LL1 - LL2 * LL2));
          A := Pi - (B + C);
          S := (LL1 - LL2) / (LL1 + LL2);
          OBSCUR := (S * S * A + B - S * Sin(C)) / Pi;
        end;
      ANNULAR:
        begin
          S := (LL1 - LL2) / (LL1 + LL2);
          OBSCUR := S * S;
        end;
      TOTAL:
        OBSCUR := 1.0;
    end;

  end;

  (* ----------------------------------------------------------------------- *)
  (* Use the Pegasus method to locate the times of 2nd and 3rd contact *)
  (* starting at the time of maximum eclipse. Search interval *)
  (* [T_MAX-DT,T_MAX] for 2nd contact, [T_MAX,T_MAX+DT] for 3rd contact. *)
  (* ----------------------------------------------------------------------- *)

  if (PHASE > PARTIAL) then
  begin

    DT := 0.25 / CENT; (* 15 min. in julian centuries *)

    Pegasus(UMBRAS_DIST, T_MAX - DT, T_MAX, EPS, OK, CONTACT_TIMES[2]);
    Pegasus(UMBRAS_DIST, T_MAX, T_MAX + DT, EPS, OK, CONTACT_TIMES[3]);

  end;

end; (* CONTACTS *)

(* --------------------------------------------------------------------------- *)
(* *)
(* POS_ANGLES *)
(* *)
(* Computes the position angles w.r.t. North (standard definition) and *)
(* with respect to the local vertical. *)
(* *)
(* ETDIFUT         Difference ephemeris time - universal time (s) *)
(* LAMBDA,RCOSPHI, *)
(* RSINPHI         Geocentric coordinates of the observer *)
(* POLYEPH         Chebychev polynomials for solar and lunar coordinates *)
(* CONTACT_TIMES   Contact times (in julian centuries since J2000 UT) *)
(* PHASE           Maximum phase of the eclipse *)
(* P               Position angle measured from north (in deg) *)
(* V               Position angle with respect to the vertical (in deg) *)
(* *)
(* --------------------------------------------------------------------------- *)

procedure POS_ANGLES(ETDIFUT: Double; LAMBDA, RCOSPHI, RSINPHI: Double; POLYEPH: TPOLYEPH;
  CONTACT_TIMES: REAL4; PHASE: PHASE_TYPE; var P, V: REAL4);

var
  INDX: integer;
  T_UT: Double;
  XI, ETA, ZETA: Double;
  XS, YS: Double;
  F1, L1, F2, L2, LL: Double;
  IJK: Mat3X;

begin

  for INDX := 1 to 4 do
  begin

    if (PHASE > PARTIAL) or ((PHASE = PARTIAL) and ((INDX = 1) or (INDX = 4))) then

    begin

      T_UT := CONTACT_TIMES[INDX];
      Bessel(T_UT, ETDIFUT, POLYEPH, IJK, XS, YS, F1, L1, F2, L2);
      Observer(T_UT, LAMBDA, RCOSPHI, RSINPHI, IJK, XI, ETA, ZETA);

      case INDX of
        1, 4:
          LL := L1 - ZETA * TN(F1);
        2, 3:
          LL := L2 - ZETA * TN(F2);
      end;

      P[INDX] := MOD360(ATN2((XS - XI) / LL, (YS - ETA) / LL));
      V[INDX] := MOD360(P[INDX] - ATN2(XI, ETA));

    end

    else

    begin
      P[INDX] := 0.0;
      V[INDX] := 0.0;
    end;

  end;

end; (* POS_ANGLES *)

(* --------------------------------------------------------------------------- *)
(* *)
(* Print out *)
(* *)
(* --------------------------------------------------------------------------- *)

procedure PRINT_RESULTS(CONTACT_TIMES, POSANG_P, POSANG_V: REAL4; T_MAX, MAG, MAX_OBSC: Double;
  MAX_PHASE: PHASE_TYPE);

var
  INDX: integer;
  DAY, MONTH, YEAR, HRS, MINS: integer;
  MJD_UT, HOUR, SECS: Double;

begin

  (* Eclipse type, magnitude and degree of ocultation *)

  writeln;

  case (MAX_PHASE) of
    NO_ECLIPSE:
      write(' The eclipse is not visible at the', ' selected site.');
    PARTIAL:
      write(' Partial eclipse', ' with M=', MAG:5:3, ' (', MAX_OBSC:4:2, ').');
    ANNULAR:
      write(' Annular eclipse', ' with M=', MAG:5:3, ' (', MAX_OBSC:4:2, ').');
    TOTAL:
      write(' Total eclipse', ' with M=', MAG:5:3, '.');
  end;

  (* Print contact times and position angles, if the eclipse *)
  (* is at least partial at the place of the observer *)

  if (MAX_PHASE > NO_ECLIPSE) then
  begin

    MJD_UT := T_MAX * 36525.0 + 51544.5;
    CalDat(MJD_UT, DAY, MONTH, YEAR, HOUR);
    DMS(HOUR, HRS, MINS, SECS);

    writeln(' Maximum at ', HRS:2, ':', MINS:2, ':', Trunc(SECS + 0.5):2, ' UT.');

    writeln;
    writeln('h  m  s [UT]    P [o]     V [o]':59);

    for INDX := 1 to 4 do
    begin
      if ((INDX = 1) or (INDX = 4) or (MAX_PHASE > PARTIAL)) then
      begin
        case INDX of
          1:
            write(' 1st contact: ');
          2:
            write(' 2nd contact: ');
          3:
            write(' 3rd contact: ');
          4:
            write(' 4th contact: ');
        end;
        MJD_UT := CONTACT_TIMES[INDX] * 36525.0 + 51544.5;
        CalDat(MJD_UT, DAY, MONTH, YEAR, HOUR);
        DMS(HOUR, HRS, MINS, SECS);
        write(YEAR:4, '/', MONTH:2, '/', DAY:2, ' ');
        write(HRS:3, MINS:3, Trunc(SECS + 0.5):3);
        writeln(POSANG_P[INDX]:15:0, POSANG_V[INDX]:10:0);
      end;
    end;

  end;

  writeln;

end; (* PRINT_RESULTS *)

(* --------------------------------------------------------------------------- *)

begin (* EclTimer main program *)

  (* Input data *)

  GetInput(TNEWMOON, ETDIFUT, LAMBDA, PHI);
  Site(PHI, RCOSPHI, RSINPHI);

  (* Chebyshev approximation of solar and lunar coordinates *)

  T_Fit_Moon(TNEWMOON - TDEVELOP, TNEWMOON + TDEVELOP, MAX_TP_DEG, POLYEPH[RAM], POLYEPH[DEM],
    POLYEPH[RM]);
  T_Fit_Sun(TNEWMOON - TDEVELOP, TNEWMOON + TDEVELOP, MAX_TP_DEG, POLYEPH[RAS], POLYEPH[DES],
    POLYEPH[RS]);

  (* Compute contact times and position angles *)

  CONTACTS(TNEWMOON, ETDIFUT, LAMBDA, RCOSPHI, RSINPHI, POLYEPH, CONTACT_TIMES, T_MAX, MAX_MAGN,
    MAX_OBSC, MAX_PHASE);

  POS_ANGLES(ETDIFUT, LAMBDA, RCOSPHI, RSINPHI, POLYEPH, CONTACT_TIMES, MAX_PHASE, POSANG_P,
    POSANG_V);

  (* Print results *)

  PRINT_RESULTS(CONTACT_TIMES, POSANG_P, POSANG_V, T_MAX, MAX_MAGN, MAX_OBSC, MAX_PHASE);

end. (* EclTimer *)
