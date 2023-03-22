unit Astro.Sphlib;

interface

uses
  Astro.Matlib,
  Astro.Pnulib;

(*-----------------------------------------------------------------------*)
(* ABERRAT: velocity vector of the Earth in equatorial coordinates       *)
(*          (in units of the velocity of light)                          *)
(*-----------------------------------------------------------------------*)
procedure ABERRAT(T: Double; var VX,VY,VZ: Double);
(*-----------------------------------------------------------------------*)
(* APPARENT  : apparent coordinates of a star                            *)
(*   PNMAT   : matrix for precession and nutation                        *)
(*   VX,VY,VZ: velocity of the earth (equatorial coord.; in units of c)  *)
(*   RA,DEC  : right ascension and declination                           *)
(*-----------------------------------------------------------------------*)
procedure APPARENT ( PNMAT: Double33; VX,VY,VZ: Double; var RA,DEC: Double);
(*-----------------------------------------------------------------------*)
(* ECLEQU: Conversion of ecliptic into equatorial coordinates            *)
(*         (T: equinox in Julian centuries since J2000)                  *)
(*-----------------------------------------------------------------------*)
procedure ECLEQU(T: Double; var X, Y, Z: Double);
(*-----------------------------------------------------------------------*)
(* ECLEQU: Conversion of ecliptic into equatorial coordinates            *)
(*         (T: equinox in Julian centuries since J2000)                  *)
(*-----------------------------------------------------------------------*)
procedure EQUECL(T: double;var X,Y,Z: double);
(* ----------------------------------------------------------------------- *)
(* EQUHOR: conversion of equatorial into horizontal coordinates *)
(* DEC  : declination (-90 deg .. +90 deg) *)
(* TAU  : hour angle (0 deg .. 360 deg) *)
(* PHI  : geographical latitude (in deg) *)
(* H    : altitude (in deg) *)
(* AZ   : azimuth (0 deg .. 360 deg, counted S->W->N->E->S) *)
(* ----------------------------------------------------------------------- *)
procedure EQUHOR(DEC, TAU, PHI: Double; var H, AZ: Double);
(*-----------------------------------------------------------------------*)
(* EQUSTD: transformation of equatorial coordinates into                 *)
(*         standard coordinates                                          *)
(*   RA0,DEC0: right ascension and declination of the optical axis (deg) *)
(*   RA,DEC:   right ascension and declination (deg)                     *)
(*   XX,YY:    standard coordinates                                      *)
(*-----------------------------------------------------------------------*)
procedure EQUSTD(RA0, DEC0, RA, DEC: Double; var XX, YY: Double);
(*-----------------------------------------------------------------------*)
(* GAUSVEC: calculation of the Gaussian vectors (P,Q,R) from             *)
(*          ecliptic orbital elements:                                   *)
(*          LAN = longitude of the ascending node                        *)
(*          INC = inclination                                            *)
(*          AOP = argument of perihelion                                 *)
(*-----------------------------------------------------------------------*)
procedure GAUSVEC(LAN, INC, AOP: Double; var PQR: Double33);
(*-----------------------------------------------------------------------*)
(* HOREQU: conversion of horizontal to equatorial coordinates            *)
(*   H,AZ : azimuth and altitude (in deg)                                *)
(*   PHI  : geographical latitude (in deg)                               *)
(*   DEC  : declination (-90 deg .. +90 deg)                             *)
(*   TAU  : hour angle (0 deg .. 360 deg)                                *)
(*-----------------------------------------------------------------------*)
procedure HOREQU(H, AZ, PHI: Double; var DEC, TAU: Double);
(*-----------------------------------------------------------------------*)
(* ORBECL: transformation of coordinates XX,YY referred to the           *)
(*         orbital plane into ecliptic coordinates X,Y,Z using           *)
(*         Gaussian vectors PQR                                          *)
(*-----------------------------------------------------------------------*)
procedure ORBECL(XX, YY: Double; PQR: Double33; var X, Y, Z: Double);
(*-----------------------------------------------------------------------*)
(* SITE:  calculates geocentric from geographic coordinates              *)
(*        RCPHI:  r * cos(phi') (geocentric; in earth radii)             *)
(*        RSPHI:  r * sin(phi') (geocentric; in earth radii)             *)
(*        PHI:    geographic latitude (deg)                              *)
(*-----------------------------------------------------------------------*)
procedure SITE ( PHI: Double; var RCPHI,RSPHI: Double);
(*-----------------------------------------------------------------------*)
(* STDEQU: transformation from standard coordinates into                 *)
(*         equatorial coordinates                                        *)
(*   RA0,DEC0: right ascension and declination of the optical axis (deg) *)
(*   XX,YY:    standard coordinates                                      *)
(*   RA,DEC:   right ascension and declination (deg)                     *)
(*-----------------------------------------------------------------------*)
procedure STDEQU ( RA0,DEC0,XX,YY: Double; var RA,DEC: Double);


implementation

(*-----------------------------------------------------------------------*)

procedure ABERRAT(T: Double; var VX, VY, VZ: Double);
const
  P2 = 6.283185307;
var
  L, CL: Double;

  function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1;
    Frac := X
  end;

begin
  L := P2 * Frac(0.27908 + 100.00214 * T);
  CL := Cos(L);
  VX := -0.994E-4 * Sin(L);
  VY := +0.912E-4 * CL;
  VZ := +0.395E-4 * CL;
end;

(*-----------------------------------------------------------------------*)

procedure APPARENT(PNMAT: Double33; VX, VY, VZ: Double; var RA, DEC: Double);
var
  X, Y, Z, R: Double;
begin
  CART(1.0, DEC, RA, X, Y, Z); (* cartesian coordinates of the star *)
  PRECART(PNMAT, X, Y, Z); (* correct for precession and nutation *)
  X := X + VX;
  Y := Y + VY;
  Z := Z + VZ; (* aberration *)
  POLAR(X, Y, Z, R, DEC, RA); (* apparent right ascension,declination *)
end;

(*-----------------------------------------------------------------------*)

procedure ECLEQU(T: Double; var X, Y, Z: Double);
var
  EPS, C, S, V: Double;
begin
  EPS := 23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0;
  C := CS(EPS);
  S := SN(EPS);
  V := +C * Y - S * Z;
  Z := +S * Y + C * Z;
  Y := V;
end;

(*-----------------------------------------------------------------------*)

procedure EQUECL(T: Double; var X, Y, Z: Double);
var
  EPS, C, S, V: Double;
begin
  EPS := 23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0;
  C := CS(EPS);
  S := SN(EPS);
  V := +C * Y + S * Z;
  Z := -S * Y + C * Z;
  Y := V;
end;

(* ----------------------------------------------------------------------- *)

procedure EQUHOR(DEC, TAU, PHI: Double; var H, AZ: Double);
var
  CS_PHI, SN_PHI, CS_DEC, SN_DEC, CS_TAU, X, Y, Z, DUMMY: Double;
begin
  CS_PHI := CS(PHI);
  SN_PHI := SN(PHI);
  CS_DEC := CS(DEC);
  SN_DEC := SN(DEC);
  CS_TAU := CS(TAU);
  X := CS_DEC * SN_PHI * CS_TAU - SN_DEC * CS_PHI;
  Y := CS_DEC * SN(TAU);
  Z := CS_DEC * CS_PHI * CS_TAU + SN_DEC * SN_PHI;
  POLAR(X, Y, Z, DUMMY, H, AZ)
end;

(*-----------------------------------------------------------------------*)

procedure EQUSTD(RA0, DEC0, RA, DEC: Double; var XX, YY: Double);
var
  C: Double;
begin
  C := CS(DEC0) * CS(DEC) * CS(RA - RA0) + SN(DEC0) * SN(DEC);
  XX := -(CS(DEC) * SN(RA - RA0)) / C;
  YY := -(SN(DEC0) * CS(DEC) * CS(RA - RA0) - CS(DEC0) * SN(DEC)) / C;
end;

(*-----------------------------------------------------------------------*)

procedure GAUSVEC(LAN, INC, AOP: Double; var PQR: Double33);
var
  C1, S1, C2, S2, C3, S3: Double;
begin
  C1 := CS(AOP);
  C2 := CS(INC);
  C3 := CS(LAN);
  S1 := SN(AOP);
  S2 := SN(INC);
  S3 := SN(LAN);
  PQR[1, 1] := +C1 * C3 - S1 * C2 * S3;
  PQR[1, 2] := -S1 * C3 - C1 * C2 * S3;
  PQR[1, 3] := +S2 * S3;
  PQR[2, 1] := +C1 * S3 + S1 * C2 * C3;
  PQR[2, 2] := -S1 * S3 + C1 * C2 * C3;
  PQR[2, 3] := -S2 * C3;
  PQR[3, 1] := +S1 * S2;
  PQR[3, 2] := +C1 * S2;
  PQR[3, 3] := +C2;
end;

(*-----------------------------------------------------------------------*)

procedure HOREQU(H, AZ, PHI: Double; var DEC, TAU: Double);
var
  CS_PHI, SN_PHI, CS_H, SN_H, CS_AZ, X, Y, Z, DUMMY: Double;
begin
  CS_PHI := CS(PHI);
  SN_PHI := SN(PHI);
  CS_H := CS(H);
  SN_H := SN(H);
  CS_AZ := CS(AZ);
  X := CS_H * SN_PHI * CS_AZ + SN_H * CS_PHI;
  Y := CS_H * SN(AZ);
  Z := SN_H * SN_PHI - CS_H * CS_PHI * CS_AZ;
  POLAR(X, Y, Z, DUMMY, DEC, TAU)
end;

(*-----------------------------------------------------------------------*)

procedure ORBECL(XX, YY: Double; PQR: Double33; var X, Y, Z: Double);
begin
  X := PQR[1, 1] * XX + PQR[1, 2] * YY;
  Y := PQR[2, 1] * XX + PQR[2, 2] * YY;
  Z := PQR[3, 1] * XX + PQR[3, 2] * YY;
end;

(*-----------------------------------------------------------------------*)

procedure SITE(PHI: Double; var RCPHI, RSPHI: Double);
const
  E2 = 0.006694; (* e**2 = f(2-f) for flattening f = 1/298.257 *)
var
  N, SNPHI: Double;
begin
  SNPHI := SN(PHI);
  N := 1.0 / Sqrt(1.0 - E2 * SNPHI * SNPHI);
  RCPHI := N * CS(PHI);
  RSPHI := (1.0 - E2) * N * SNPHI;
end;

(*-----------------------------------------------------------------------*)

procedure STDEQU(RA0, DEC0, XX, YY: Double; var RA, DEC: Double);
begin
  RA := RA0 + ATN(-XX / (CS(DEC0) - YY * SN(DEC0)));
  DEC := ASN((SN(DEC0) + YY * CS(DEC0)) / Sqrt(1.0 + XX * XX + YY * YY));
end;



end.
