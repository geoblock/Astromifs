unit Apc.Kepler;

interface

uses
  System.Math,
  Apc.Mathem,
  Apc.Spheric;

(*-----------------------------------------------------------------------*)
(* EccAnom: calculation of the eccentric anomaly E=EccAnom(MAN,ECC)      *)
(*          from the mean anomaly MAN and the eccentricity ECC.          *)
(*          (solution of Kepler's equation by Newton's method)           *)
(*          (E, MAN in degrees)                                          *)
(*-----------------------------------------------------------------------*)
function EccAnom(MAN, ECC: Double): Double;
(*-----------------------------------------------------------------------*)
(* Element: calculates orbital elements from two positions               *)
(*          for elliptic and hyperbolic orbits                           *)
(*                                                                       *)
(*   JDA,JDB: time of passage of points A and B (Julian Date)            *)
(*   Ra, RB : position vectors of points A and B                         *)
(*   TP     : perihelion time (in Julian centuries since J2000)          *)
(*   Q      : perihelion distance                                        *)
(*   ECC    : eccentricity                                               *)
(*   INC    : inclination (in deg)                                       *)
(*   LAN    : longitude of the ascending node (in deg)                   *)
(*   AOP    : argument of perihelion (in deg)                            *)
(*-----------------------------------------------------------------------*)
procedure Element(JDA, JDB: Double; Ra, RB: Vector; var TP, Q, ECC, INC, LAN, AOP: Double);
(*-----------------------------------------------------------------------*)
(* Ellip: calculation of position and velocity vector                    *)
(*        for elliptic orbits                                            *)
(*                                                                       *)
(*        M    mean anomaly in degrees       X,Y    position   (in AU)   *)
(*        A    semimajor axis                VX,VY  velocity (in AU/day) *)
(*        ECC  eccentricity                                              *)
(*-----------------------------------------------------------------------*)
procedure Ellip(M, A, ECC: Double; var X, Y, VX, VY: Double);
(*-----------------------------------------------------------------------*)
(* Find_ETA: determines the sector/triangle ratio                        *)
(*           from two positions and the time difference                  *)
(*-----------------------------------------------------------------------*)
function Find_ETA(Ra, RB: Vector; TAU: Double): Double;
(*-----------------------------------------------------------------------*)
(* HypAnom: calculation of the eccentric anomaly H=HypAnom(MH,ECC) from  *)
(*          mean anomaly MH and eccentricity ECC for                     *)
(*          hyperbolic orbits                                            *)
(*-----------------------------------------------------------------------*)
function HypAnom(MH, ECC: Double): Double;
(*-----------------------------------------------------------------------*)
(* Hyperb: calculation of the position and velocity vector               *)
(*         for hyperbolic orbits                                         *)
(*                                                                       *)
(*   T0   time of perihelion passage             X,Y   position          *)
(*   T    time                                   VX,VY velocity          *)
(*   A    semimajor axis (arbitrary sign)                                *)
(*   ECC  eccentricity                                                   *)
(*   (T0,T in julian centuries since J2000)                              *)
(*-----------------------------------------------------------------------*)
procedure Hyperb(T0, T, A, ECC: Double; var X, Y, VX, VY: Double);
(*-----------------------------------------------------------------------*)
(* Kepler: calculation of position and velocity for unperturbed          *)
(*         elliptic, parabolic and hyperbolic orbits                     *)
(*                                                                       *)
(*        T0   time of perihelion passage         X,Y,Z     position     *)
(*        T    time                               VX,VY,VZ  velocity     *)
(*        Q    perihelion distance                                       *)
(*        ECC  eccentricity                                              *)
(*        PQR  matrix of Gaussian vectors                                *)
(*        (T0,T in Julian centuries since J2000)                         *)
(*-----------------------------------------------------------------------*)
procedure Kepler(T0, T, Q, ECC: Double; PQR: Double33; var X, Y, Z, VX, VY, VZ: Double);
(*-----------------------------------------------------------------------*)
(* Parab: calculation of position and velocity for                       *)
(*        parabolic and near parabolic orbits according to Stumpff       *)
(*                                                                       *)
(*        T0   time of perihelion passage         X,Y    position        *)
(*        T    time                              VX,VY  velocity         *)
(*        Q    perihelion distance                                       *)
(*        ECC  eccentricity                                              *)
(*        (T0,T in julian centuries since J2000)                         *)
(*-----------------------------------------------------------------------*)
procedure Parab(T0, T, Q, ECC: Double; var X, Y, VX, VY: Double);
(*-----------------------------------------------------------------------*)
(*  Stumpff: calculation of Stumpff's functions C1 = sin(E)/E ,          *)
(*           C2 = (1-cos(E))/(E**2) and C3 = (E-sin(E))/(E**3)           *)
(*           for argument E2=E**2                                        *)
(*           (E: eccentric anomaly in radian)                            *)
(*-----------------------------------------------------------------------*)
procedure Stumpff(E2: Double; var C1, C2, C3: Double);
(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* XYZKep: conversion of the state vector into Keplerian elements        *)
(*         for elliptical orbits                                         *)
(*                                                                       *)
(*  X,Y,Z    : Position [AU]                                             *)
(*  VX,VY,VZ : Velocity [AU/d]                                           *)
(*  AX       : Semi-major axis [AU]                                      *)
(*  ECC      : Eccentricity                                              *)
(*  INC      : Inclination [deg]                                         *)
(*  LAN      : Longitude of the ascending node [deg]                     *)
(*  AOP      : Argument of perihelion [deg]                              *)
(*  MA       : Mean anomaly [deg]                                        *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure XYZKep(X, Y, Z, VX, VY, VZ: Double; var AX, ECC, INC, LAN, AOP, MA: Double);


//========================================================================
implementation
//========================================================================

(*-----------------------------------------------------------------------*)
function EccAnom(MAN, ECC: Double): Double;
const
  TWOPI = 6.283185308;
  RAD = 0.0174532925199433;
  EPS = 1E-11;
  MAXIT = 15;
var
  M, E, F: Double;
  I: integer;
begin
  M := MAN / 360.0;
  M := TWOPI * (M - Trunc(M));
  if M < 0 then
    M := M + TWOPI;
  if (ECC < 0.8) then
    E := M
  else
    E := Pi;
  F := E - ECC * Sin(E) - M;
  I := 0;
  while((abs(F) > EPS) and (I < MAXIT))
  do
  begin
    E := E - F / (1.0 - ECC * Cos(E));
    F := E - ECC * Sin(E) - M;
    I := I + 1;
  end;
  EccAnom := E / RAD;
  if (I = MAXIT) then
    writeln(' convergence problems in EccAnom');
end;

(*-----------------------------------------------------------------------*)

procedure Element(JDA, JDB: Double; Ra, RB: Vector; var TP, Q, ECC, INC, LAN, AOP: Double);
const
  KGAUSS = 0.01720209895;
  RAD = 0.01745329252; (* 180/pi *)

var
  TAU, ETA, P, AX, N, NY, E, M, U: Double;
  SA, SB, S0, FAC, DUMMY, SHH: Double;
  COS_DNY, SIN_DNY, ECOS_NY, ESIN_NY: Double;
  EA, R0, E0, R: Vector;
  I: Index;

begin
  (* calculate vector R0 (fraction of RB perpendicular to Ra) *)
  (* and the magnitudes of Ra, RB and R0 *)
  SA := Norm(Ra);
  for I := X to Z do
    EA[I] := Ra[I] / SA;
  SB := Norm(RB);
  FAC := Dot(RB, EA);
  for I := X to Z do
    R0[I] := RB[I] - FAC * EA[I];
  S0 := Norm(R0);
  for I := X to Z do
    E0[I] := R0[I] / S0;

  (* inclination and ascending node *)
  Cross(EA, E0, R);
  Polar(-R[Y], R[X], R[Z], DUMMY, INC, LAN);
  INC := 90.0 - INC;
  U := ATN2((+E0[X] * R[Y] - E0[Y] * R[X]), (-EA[X] * R[Y] + EA[Y] * R[X]));
  if INC = 0.0 then
    U := ATN2(Ra[Y], Ra[X]);

  (* semilatus rectum p *)
  TAU := KGAUSS * abs(JDB - JDA);
  ETA := Find_ETA(Ra, RB, TAU);
  P := SA * S0 * ETA / TAU;
  P := P * P;

  (* eccentricity, true anomaly and longitude of perihelion *)
  COS_DNY := FAC / SB;
  SIN_DNY := S0 / SB;
  ECOS_NY := P / SA - 1.0;
  ESIN_NY := (ECOS_NY * COS_DNY - (P / SB - 1.0)) / SIN_DNY;
  Polar(ECOS_NY, ESIN_NY, 0.0, ECC, DUMMY, NY);
  AOP := U - NY;
  while (AOP < 0.0) do
    AOP := AOP + 360.0;

  (* perihelion distance, semimajor axis and mean daily motion *)
  Q := P / (1.0 + ECC);
  AX := Q / (1.0 - ECC);
  N := KGAUSS / Sqrt(abs(AX * AX * AX));

  (* mean anomaly and time of perihelion passage *)
  if (ECC < 1.0) then
  begin
    E := ATN2(Sqrt((1.0 - ECC) * (1.0 + ECC)) * ESIN_NY, ECOS_NY + ECC * ECC);
    E := RAD * E;
    M := E - ECC * Sin(E);;
  end
  else
  begin
    SHH := Sqrt((ECC - 1.0) * (ECC + 1.0)) * ESIN_NY / (ECC + ECC * ECOS_NY);
    M := ECC * SHH - LN(SHH + Sqrt(1.0 + SHH * SHH))
  end;
  TP := ((JDA - M / N) - 2451545.0) / 36525.0;
end;

(* ----------------------------------------------------------------------- *)

procedure Ellip(M, A, ECC: Double; var X, Y, VX, VY: Double);
const
  KGAUSS = 0.01720209895;
var
  K, E, C, S, FAC, RHO: Double;
begin
  K := KGAUSS / Sqrt(A);
  E := EccAnom(M, ECC);
  C := CS(E);
  S := SN(E);
  FAC := Sqrt((1.0 - ECC) * (1 + ECC));
  RHO := 1.0 - ECC * C;
  X := A * (C - ECC);
  Y := A * FAC * S;
  VX := -K * S / RHO;
  VY := K * FAC * C / RHO;
end;

(* ----------------------------------------------------------------------- *)

function Find_ETA(Ra, RB: Vector; TAU: Double): Double;

const
  DELTA = 1.0E-9;
  MAXIT = 30;

var
  KAPPA, M, L, SA, SB, ETA_MIN, ETA1, ETA2, F1, F2, D_ETA: Double;
  I: integer;

  (* F(eta) = 1 - eta + (m/eta**2)*W(m/eta**2-l) *)
  (*sub*)function F(ETA, M, L: Double): Double;
  const
    EPS = 1.0E-10;
  var
    W, WW, A, S, N, G, E: Double;
  begin
    W := M / (ETA * ETA) - L;
    if (abs(W) < 0.1) then (* series expansion *)
    begin
      A := 4.0 / 3.0;
      WW := A;
      N := 0.0;
      repeat
        N := N + 1;
        A := A * W * (N + 2.0) / (N + 1.5);
        WW := WW + A;
      until abs(A) < EPS;
    end
    else if (W > 0) then (* W=(2g-sin2g)/(sin(g)**3), g=2*arcsin(sqrt(w)) *)
    begin
      G := 2.0 * ArcTan(Sqrt(W / (1.0 - W)));
      S := Sin(G);
      WW := (2.0 * G - Sin(2.0 * G)) / (S * S * S);
    end
    else (* W=(sinh2g-2g)/(sinh(g)**3), g=2*arsinh(sqrt(-w)) *)
    begin
      G := 2.0 * LN(Sqrt(-W) + Sqrt(1.0 - W));
      E := EXP(G);
      S := 0.5 * (E - 1.0 / E);
      E := E * E;
      WW := (0.5 * (E - 1.0 / E) - 2.0 * G) / (S * S * S);
    end;
    F := 1.0 - ETA + (W + L) * WW;
  end; (* Find_ETA.F *)

begin
  SA := Norm(Ra);
  SB := Norm(RB);
  KAPPA := Sqrt(2.0 * (SA * SB + Dot(Ra, RB)));
  M := TAU * TAU / (KAPPA * KAPPA * KAPPA);
  L := (SA + SB) / (2.0 * KAPPA) - 0.5;
  ETA_MIN := Sqrt(M / (L + 1.0));

  (* start with Hansen's approximation *)
  ETA2 := (12.0 + 10.0 * Sqrt(1.0 + (44.0 / 9.0) * M / (L + 5.0 / 6.0))) / 22.0;
  ETA1 := ETA2 + 0.1;
  F1 := F(ETA1, M, L);
  F2 := F(ETA2, M, L);
  I := 0;

  (* secant method *)
  while ((abs(F2 - F1) > DELTA) and (I < MAXIT)) do
  begin
    D_ETA := -F2 * (ETA2 - ETA1) / (F2 - F1);
    ETA1 := ETA2;
    F1 := F2;
    while (ETA2 + D_ETA <= ETA_MIN) do
      D_ETA := 0.5 * D_ETA;
    ETA2 := ETA2 + D_ETA;
    F2 := F(ETA2, M, L);
    I := I + 1;
  end;
  if (I = MAXIT) then
    writeln(' convergence problems in Find_ETA');
  Find_ETA := ETA2;

end;

(* ----------------------------------------------------------------------- *)

function HypAnom(MH, ECC: Double): Double;
const
  EPS = 1E-10;
  MAXIT = 15;
var
  H, F, EXH, SINHH, COSHH: Double;
  I: integer;
begin
  H := LN(2.0 * abs(MH) / ECC + 1.8);
  if (MH < 0.0) then
    H := -H;
  EXH := EXP(H);
  SINHH := 0.5 * (EXH - 1.0 / EXH);
  COSHH := 0.5 * (EXH + 1.0 / EXH);
  F := ECC * SINHH - H - MH;
  I := 0;
  while ((abs(F) > EPS * (1.0 + abs(H + MH))) and (I < MAXIT)) do
  begin
    H := H - F / (ECC * COSHH - 1.0);
    EXH := EXP(H);
    SINHH := 0.5 * (EXH - 1.0 / EXH);
    COSHH := 0.5 * (EXH + 1.0 / EXH);
    F := ECC * SINHH - H - MH;
    I := I + 1;
  end;
  HypAnom := H;
  if (I = MAXIT) then
    writeln(' convergence problems in HypAnom');
end;

(*-----------------------------------------------------------------------*)

procedure Hyperb(T0, T, A, ECC: Double; var X, Y, VX, VY: Double);
const
  KGAUSS = 0.01720209895;
var
  K, MH, H, EXH, COSHH, SINHH, RHO, FAC: Double;
begin
  A := abs(A);
  K := KGAUSS / Sqrt(A);
  MH := K * 36525.0 * (T - T0) / A;
  H := HypAnom(MH, ECC);
  EXH := EXP(H);
  COSHH := 0.5 * (EXH + 1.0 / EXH);
  SINHH := 0.5 * (EXH - 1.0 / EXH);
  FAC := Sqrt((ECC + 1.0) * (ECC - 1.0));
  RHO := ECC * COSHH - 1.0;
  X := A * (ECC - COSHH);
  Y := A * FAC * SINHH;
  VX := -K * SINHH / RHO;
  VY := K * FAC * COSHH / RHO;
end;

(* ----------------------------------------------------------------------- *)

procedure Kepler(T0, T, Q, ECC: Double; PQR: Double33; var X, Y, Z, VX, VY, VZ: Double);
const
  M0 = 5.0;
  EPS = 0.1;
  KGAUSS = 0.01720209895;
  DEG = 57.29577951;
var
  M, DELTA, TAU, INVAX, XX, YY, VVX, VVY: Double;
begin
  DELTA := abs(1.0 - ECC);
  INVAX := DELTA / Q;
  TAU := KGAUSS * 36525.0 * (T - T0);
  M := DEG * TAU * Sqrt(INVAX * INVAX * INVAX);
  if ((M < M0) and (DELTA < EPS)) then
    Parab(T0, T, Q, ECC, XX, YY, VVX, VVY)
  else if (ECC < 1.0) then
    Ellip(M, 1.0 / INVAX, ECC, XX, YY, VVX, VVY)
  else
    Hyperb(T0, T, 1.0 / INVAX, ECC, XX, YY, VVX, VVY);
  Orb2Ecl(XX, YY, PQR, X, Y, Z);
  Orb2Ecl(VVX, VVY, PQR, VX, VY, VZ);
end;

(* ----------------------------------------------------------------------- *)

procedure Parab(T0, T, Q, ECC: Double; var X, Y, VX, VY: Double);
const
  EPS = 1E-9;
  KGAUSS = 0.01720209895;
  MAXIT = 15;
var
  E2, E20, FAC, C1, C2, C3, K, TAU, A, U, U2: Double;
  R: Double;
  I: integer;
begin
  E2 := 0.0;
  FAC := 0.5 * ECC;
  I := 0;
  K := KGAUSS / Sqrt(Q * (1.0 + ECC));
  TAU := KGAUSS * 36525.0 * (T - T0);
  repeat
    I := I + 1;
    E20 := E2;
    A := 1.5 * Sqrt(FAC / (Q * Q * Q)) * TAU;
    A := Cubr(Sqrt(A * A + 1.0) + A);
    U := A - 1.0 / A;
    U2 := U * U;
    E2 := U2 * (1.0 - ECC) / FAC;
    Stumpff(E2, C1, C2, C3);
    FAC := 3.0 * ECC * C3;
  until (abs(E2 - E20) < EPS) or (I > MAXIT);
  if (I = MAXIT) then
    writeln(' convergence problems in Parab');
  R := Q * (1.0 + U2 * C2 * ECC / FAC);
  X := Q * (1.0 - U2 * C2 / FAC);
  VY := K * (X / R + ECC);
  Y := Q * Sqrt((1.0 + ECC) / FAC) * U * C1;
  VX := -K * Y / R;
end;

(*-----------------------------------------------------------------------*)

procedure Stumpff(E2: Double; var C1, C2, C3: Double);
const
  EPS = 1E-12;
var
  N, ADD: Double;
begin
  C1 := 0.0;
  C2 := 0.0;
  C3 := 0.0;
  ADD := 1.0;
  N := 1.0;
  repeat
    C1 := C1 + ADD;
    ADD := ADD / (2.0 * N);
    C2 := C2 + ADD;
    ADD := ADD / (2.0 * N + 1.0);
    C3 := C3 + ADD;
    ADD := -E2 * ADD;
    N := N + 1.0;
  until abs(ADD) < EPS;
end;

(*-----------------------------------------------------------------------*)

procedure XYZKep(X, Y, Z, VX, VY, VZ: Double; var AX, ECC, INC, LAN, AOP, MA: Double);
const
  DEG = 57.29577951308; (* Conversion from radian to degrees *)
  KGAUSS = 0.01720209895; (* Gaussian gravitational constant *)
var
  HX, HY, HZ, H, R, V2: Double;
  GM, C, S, E2, EA, U, NU: Double;
begin
  HX := Y * VZ - Z * VY; (* Areal velocity *)
  HY := Z * VX - X * VZ;
  HZ := X * VY - Y * VX;
  H := Sqrt(HX * HX + HY * HY + HZ * HZ);

  LAN := ATN2(HX, -HY); (* Long. ascend. node *)
  INC := ATN2(Sqrt(HX * HX + HY * HY), HZ); (* Inclination *)
  U := ATN2(Z * H, -X * HY + Y * HX); (* Arg. of latitude *)

  GM := KGAUSS * KGAUSS;
  R := Sqrt(X * X + Y * Y + Z * Z); (* Distance *)
  V2 := VX * VX + VY * VY + VZ * VZ; (* Velocity squared *)

  AX := 1.0 / (2.0 / R - V2 / GM); (* Semi-major axis *)

  C := 1.0 - R / AX; (* e*cos(E) *)
  S := (X * VX + Y * VY + Z * VZ) / (Sqrt(AX) * KGAUSS); (* e*sin(E) *)

  E2 := C * C + S * S;
  ECC := Sqrt(E2); (* Eccentricity *)
  EA := ATN2(S, C); (* Eccentric anomaly *)

  MA := EA - S * DEG; (* Mean anomaly *)

  NU := ATN2(Sqrt(1.0 - E2) * S, C - E2); (* True anomaly *)
  AOP := U - NU; (* Arg. of perihelion *)

  if (LAN < 0.0) then
    LAN := LAN + 360.0;
  if (AOP < 0.0) then
    AOP := AOP + 360.0;
  if (MA < 0.0) then
    MA := MA + 360.0;
end;

end.
