unit Apc.Physic;

interface

uses
  Apc.Mathem,
  Apc.Spheric,
  Apc.Planets;

type
  RotationType = (Direct, Retrograde);
  SystemType = (SYS_I, SYS_II, SYS_III);

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* PosAng: computes the position angle from given directions             *)
(*                                                                       *)
(*   X,Y,Z      Coordinates of the planet relative to the observer       *)
(*   DX,DY,DZ   Direction vector                                         *)
(*   PosAng     Position angle (0<=PosAng<360deg)                        *)
(*                                                                       *)
(* Both vectors must be given in a common coordinate system              *)
(* (e.g. mean equator and equinox of date)                               *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
function PosAng(X, Y, Z, DX, DY, DZ: Double): Double;
(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Shape:  returns the equatorial radii and flattening of planets        *)
(*                                                                       *)
(*   Planet   Name of the planet                                         *)
(*   R_EQU    Equatorial radius (km)                                     *)
(*   FL       Geometric flattening                                       *)
(*                                                                       *)
(* Source: M.E. Davies et al., Report of the IAU Working Group on        *)
(*         Cartographic Coordinates and Rotational Elements, 1982        *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure Shape(Planet: PlanetType; var R_EQU, FL: Double);
(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Orient:  returns the elements describing the planetocentric           *)
(*          coordinate system of a planet                                *)
(*                                                                       *)
(*   Planet   Name of the planet                                         *)
(*   SYSTEM   System of rotation (I, II or III)                          *)
(*   T        Time in Julian centuries since J2000  (ET or TDB/TDT)      *)
(*   A        Right ascension of the axis of rotation                    *)
(*   D        Declination of the axis of rotation                        *)
(*   W        Orientation of the prime meridian with respect to the      *)
(*            intersection of the Earth's equator of J2000 and the       *)
(*            planetary equator of date                                  *)
(*   SENSE    Sense of rotation (Direct or Retrograde)                   *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure Orient(Planet: PlanetType; SYSTEM: SystemType; T: Double; var A, D, W: Double;
  var SENSE: RotationType);
(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Rotation: computes the rotation parameters of a planet                *)
(*                                                                       *)
(*   X,Y,Z    Geocentric equatorial coordinates of the planet (AU)       *)
(*   A        Right ascension of the axis of rotation (J2000)            *)
(*   D        Declination of the axis of rotation (J2000)                *)
(*   W        Orientation of the prime meridian with respect to the      *)
(*            intersection of the Earth's equator of J2000 and the       *)
(*            planetary equator of date                                  *)
(*   SENSE    Sense of rotation (Direct or Retrograde)                   *)
(*   FLATT    Geometric flattening of the planet                         *)
(*   AX,AY,AZ Rotation axis unit vector (J2000)                          *)
(*   LONG     Planetographic longitude of the Earth (deg)                *)
(*   LAT      Planetographic latitude of the Earth (deg)                 *)
(*   Dec      Planetocentric latitude of the Earth (deg)                 *)
(*                                                                       *)
(* The heliocentric coordinates of the planet may be substituted for     *)
(* X,Y,Z to obtain the planetographic coordinates of the Sun.            *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure Rotation(X, Y, Z, A, D, W: Double; SENSE: RotationType; FLATT: Double;
  var AX, AY, AZ, LONG, LAT, Dec: Double);
(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Bright: Computes the apparent magnitude of a planet                   *)
(*                                                                       *)
(*  Planet     Name of the planet                                        *)
(*  R          Heliocentric distance of the planet (AU)                  *)
(*  DELTA      Distance of the planet from the observer (AU)             *)
(*  PHI        Phase angle (deg)                                         *)
(*  Dec        Planetocentric latitude of the observer (deg)             *)
(*  DLONG      Difference of the planetocentric longitudes of the Sun    *)
(*             and the observer (deg)                                    *)
(*                                                                       *)
(* Magnitudes V(1,0) from Astronomical Almanac 1984. Dec and DLONG are   *)
(* only required to compute the apparent brightness of Saturn, which     *)
(* depends on the ring orientation.                                      *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
function Bright(Planet: PlanetType; R, DELTA, PHI, Dec, DLONG: Double): Double;
(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Illum: Computes the illumination parameters of a planet               *)
(*                                                                       *)
(*   X,Y,Z      Heliocentric coordinates of the planet                   *)
(*   XE,YE,ZE   Heliocentric coordinates of the Earth                    *)
(*   R          Heliocentric distance of the planet                      *)
(*   D          Geocentric distance of the planet                        *)
(*   ELONG      Elongation (deg)                                         *)
(*   PHI        Phase angle (deg)                                        *)
(*   K          Phase                                                    *)
(*                                                                       *)
(* Note: All coordinates must refer to the same coordinate system.       *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure Illum(X, Y, Z, XE, YE, ZE: Double; var R, D, ELONG, PHI, K: Double);


implementation

(* ----------------------------------------------------------------------- *)
function PosAng(X, Y, Z, DX, DY, DZ: Double): Double;
var
  C, S, PHI: Double;
begin
  C := ((-X * Z) * DX + (-Y * Z) * DY + (X * X + Y * Y) * DZ) / SQRT(X * X + Y * Y + Z * Z);
  S := ((-Y) * DX + (X) * DY);
  PHI := ATN2(S, C);
  if (PHI < 0.0) then
    PosAng := PHI + 360.0
  else
    PosAng := PHI
end; // PosAng

(* ----------------------------------------------------------------------- *)
procedure Shape(Planet: PlanetType; var R_EQU, FL: Double);
begin
  case Planet of
    Mercury:
      begin
        R_EQU := 2439.0;
        FL := 0.0
      end;
    Venus:
      begin
        R_EQU := 6051.0;
        FL := 0.0
      end;
    Earth:
      begin
        R_EQU := 6378.14;
        FL := 0.00335281
      end;
    Mars:
      begin
        R_EQU := 3393.4;
        FL := 0.0051865
      end;
    Jupiter:
      begin
        R_EQU := 71398.0;
        FL := 0.0648088
      end;
    Saturn:
      begin
        R_EQU := 60000.0;
        FL := 0.1076209
      end;
    Uranus:
      begin
        R_EQU := 25400.0;
        FL := 0.030
      end;
    Neptune:
      begin
        R_EQU := 24300.0;
        FL := 0.0259
      end;
    Pluto:
      begin
        R_EQU := 1500.0;
        FL := 0.0
      end;
  end;
end;  // Shape

(*-----------------------------------------------------------------------*)

procedure Orient(Planet: PlanetType; SYSTEM: SystemType; T: Double; var A, D, W: Double;
  var SENSE: RotationType);
var
  TD: Double;

begin
  (* Compute right ascension and declination of the axis of *)
  (* rotation with respect to the equator and equinox of J2000 *)
  case Planet of
    Mercury:
      begin
        A := 281.02 - 0.033 * T;
        D := 61.45 - 0.005 * T
      end;
    Venus:
      begin
        A := 272.78;
        D := 67.21
      end;
    Earth:
      begin
        A := 0.00 - 0.641 * T;
        D := 90.00 - 0.557 * T
      end;
    Mars:
      begin
        A := 317.681 - 0.108 * T;
        D := 52.886 - 0.061 * T
      end;
    Jupiter:
      begin
        A := 268.05 - 0.009 * T;
        D := 64.49 + 0.003 * T
      end;
    Saturn:
      begin
        A := 40.66 - 0.036 * T;
        D := 83.52 - 0.004 * T
      end;
    Uranus:
      begin
        A := 257.43;
        D := -15.10
      end;
    Neptune:
      begin
        A := 295.33;
        D := 40.65
      end;
    Pluto:
      begin
        A := 311.63;
        D := 4.18
      end
  end;

  (* Compute orientation of the prime meridian *)

  TD := 36525.0 * T;
  case Planet of
    Mercury:
      W := 329.710 + 6.1385025 * TD;
    Venus:
      W := 159.910 - 1.4814205 * TD;
    Earth:
      W := 100.21 + 360.9856123 * TD;
    Mars:
      W := 176.655 + 350.8919830 * TD;
    Jupiter:
      case SYSTEM of
        SYS_I:
          W := 67.10 + 877.900 * TD;
        SYS_II:
          W := 43.30 + 870.270 * TD;
        SYS_III:
          W := 284.95 + 870.536 * TD
      end;
    Saturn:
      case SYSTEM of
        SYS_I, SYS_II:
          W := 227.2037 + 844.3000000 * TD;
        SYS_III:
          W := 38.90 + 810.7939024 * TD
      end;
    Uranus:
      W := 261.620 - 554.913 * TD; (* System III *)
    Neptune:
      W := 107.210 + 468.75 * TD; (* System III *)
    Pluto:
      W := 252.660 - 56.364 * TD
  end;
  W := W / 360.0;
  W := 360.0 * (W - Trunc(W));

  (* Define sense of rotation *)
  case Planet of
    Mercury, Earth, Mars, Jupiter, Saturn, Neptune:
      SENSE := Direct;
    Venus, Uranus, Pluto:
      SENSE := Retrograde;
  end;
end; // Orient

(*-----------------------------------------------------------------------*)
procedure Rotation(X, Y, Z, A, D, W: Double; SENSE: RotationType; FLATT: Double;
  var AX, AY, AZ, LONG, LAT, Dec: Double);
var
  SX, SY, SZ, T: Double;
  E: Double33;
begin

  (* Compute unit vectors E(*,1) (intersection of prime meridian and *)
  (* planetary equator), E(*,3) (parallel to the rotation axis) and *)
  (* E(*,2) perpendicular to E(*,1) and E(*,3)) *)
   GaussVec(90.0 + A, 90.0 - D, W, E);

  (* Copy rotation axis unit vector *)
   AX := E[1, 3];
  AY := E[2, 3];
  AZ := E[3, 3];

  (* Compute planetocentric latitude and longitude *)
   SX := -(E[1, 1] * X + E[2, 1] * Y + E[3, 1] * Z);
  SY := -(E[1, 2] * X + E[2, 2] * Y + E[3, 2] * Z);
  SZ := -(E[1, 3] * X + E[2, 3] * Y + E[3, 3] * Z);

  T := SZ / SQRT(SX * SX + SY * SY);

  Dec := ATN(T);
  LONG := ATN2(SY, SX);

  (* Compute planetographic latitude and longitude *)
   if (SENSE = Direct) then
    LONG := -LONG;
  if LONG < 0.0 then
    LONG := LONG + 360.0;
  LAT := ATN(T / ((1.0 - FLATT) * (1.0 - FLATT)));

end; // Rotation

(* ----------------------------------------------------------------------- *)
function Bright(Planet: PlanetType; R, DELTA, PHI, Dec, DLONG: Double): Double;
const
  LN10 = 2.302585093; (* Natural logarithm of 10 *)
var
  P, SD, DL, MAG: Double;
begin
  P := PHI / 100.0;
  case Planet of
    Mercury:
      MAG := -0.42 + (3.80 - (2.73 - 2.00 * P) * P) * P;
    Venus:
      MAG := -4.40 + (0.09 + (2.39 - 0.65 * P) * P) * P;
    Earth:
      MAG := -3.86;
    Mars:
      MAG := -1.52 + 1.6 * P;
    Jupiter:
      MAG := -9.40 + 0.5 * P;
    Saturn:
      begin
        SD := ABS(SN(Dec));
        DL := ABS(DLONG / 100.0);
        if (DL > 1.8) then
          DL := ABS(DL - 3.6);
        MAG := -8.88 - 2.60 * SD + 1.25 * SD * SD + 4.40 * DL;
      end;
    Uranus:
      MAG := -7.19;
    Neptune:
      MAG := -6.87;
    Pluto:
      MAG := -1.0;
  end;
  Bright := MAG + 5.0 * LN(R * DELTA) / LN10;
end; // Bright

(* ----------------------------------------------------------------------- *)

procedure Illum(X, Y, Z, XE, YE, ZE: Double; var R, D, ELONG, PHI, K: Double);
var
  XP, YP, ZP, RE, C_PHI: Double;
begin
  (* Compute the planet's geocentric position *)
   XP := X - XE;
  YP := Y - YE;
  ZP := Z - ZE;

  (* Compute the distances in the Sun-Earth-planet triangle *)
   R := SQRT(X * X + Y * Y + Z * Z); (* Sun-planet distance *)
  RE := SQRT(XE * XE + YE * YE + ZE * ZE); (* Sun-Earth  distance *)
  D := SQRT(XP * XP + YP * YP + ZP * ZP); (* Earth-planet distance *)

  (* Compute elongation, phase angle and phase *)
  ELONG := ACS((D * D + RE * RE - R * R) / (2.0 * D * RE));
  C_PHI := (D * D + R * R - RE * RE) / (2.0 * D * R);
  PHI := ACS(C_PHI);
  K := 0.5 * (1.0 + C_PHI);
end; // Illum

end.
