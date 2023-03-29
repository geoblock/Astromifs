unit Apc.Moon;

interface

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric;


(*-----------------------------------------------------------------------*)
(* MiniMoon: low precision lunar coordinates (approx. 5'/1')            *)
(*            T  : time in Julian centuries since J2000                  *)
(*                 ( T=(JD-2451545)/36525 )                              *)
(*            Ra : right ascension (in h; equinox of date)               *)
(*            Dec: declination (in deg; equinox of date)                 *)
(*-----------------------------------------------------------------------*)
procedure MiniMoon(T: Double; var Ra, Dec: Double);

(*-----------------------------------------------------------------------*)
(* Moon: analytical lunar theory by E.W.Brown (Improved Lunar Ephemeris) *)
(*       with an accuracy of approx. 1"                                  *)
(*                                                                       *)
(*       T:      time in Julian centuries since J2000 (Ephemeris Time)   *)
(*               (T=(JD-2451545.0)/36525.0)                              *)
(*       LAMBDA: geocentric ecliptic longitude (equinox of date)         *)
(*       BETA:   geocentric ecliptic latitude  (equinox of date)         *)
(*       R:      geocentric distance (in Earth radii)                    *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure Moon(T: Double; var LAMBDA, BETA, R: Double);

(*-----------------------------------------------------------------------*)
(* MoonEqu: geocentric equatorial coordinates of the Moon                *)
(*          referred to the true equinox of date                         *)
(*   T   time in Julian centuries ephemeris time since J2000             *)
(*       ( T = (JD-2451545.0)/36525 )                                    *)
(*   Ra  right ascension (deg)                                           *)
(*   Dec declination (deg)                                               *)
(*   R   distance (in earth radii)                                       *)
(*-----------------------------------------------------------------------*)
procedure MoonEqu(T: Double; var Ra, Dec, R: Double);

(*-----------------------------------------------------------------------*)
(* T_Fit_Moon: approximates the equatorial coordinates                   *)
(*             of the Moon by Chebyshev expansions for a                 *)
(*             given period of time of at most one month                 *)
(*                                                                       *)
(*  TA     : first date (in Julian centuries since J2000)                *)
(*  TB     : last date ( TB < TA + 1 month )                             *)
(*  N      : highest order                                               *)
(*  RA_POLY: coefficients for right ascension                            *)
(*  DE_POLY: coefficients for declination                                *)
(*  R_POLY : coefficients for geocentric distance                        *)
(*-----------------------------------------------------------------------*)
procedure T_Fit_Moon(TA, TB: Double; N: integer; var RA_POLY, DE_POLY, R_POLY: TPolynomCheb);

//------------------------------------------------------------------------
implementation
//------------------------------------------------------------------------

procedure MiniMoon(T: Double; var Ra, Dec: Double);
const
  P2 = 6.283185307;
  ARC = 206264.8062;
  COSEPS = 0.91748;
  SINEPS = 0.39778; (* cos/sin(obliquity ecliptic) *)
var
  L0, L, LS, F, D, H, S, N, DL, CB: Double;
  L_MOON, B_MOON, V, W, X, Y, Z, RHO: Double;
  function Frac(X: Double): Double;
  (* with some compilers it may be necessary to replace *)
  (* TRUNC by LongTrunc oder INT if T<-24! *)
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1;
    Frac := X
  end;

begin
  (* mean elements of lunar orbit *)
  L0 := Frac(0.606433 + 1336.855225 * T); (* mean longitude Moon (in rev) *)
  L := P2 * Frac(0.374897 + 1325.552410 * T); (* mean anomaly of the Moon *)
  LS := P2 * Frac(0.993133 + 99.997361 * T); (* mean anomaly of the Sun *)
  D := P2 * Frac(0.827361 + 1236.853086 * T); (* diff. longitude Moon-Sun *)
  F := P2 * Frac(0.259086 + 1342.227825 * T); (* mean argument of latitude *)
  DL := +22640 * Sin(L) - 4586 * Sin(L - 2 * D) + 2370 * Sin(2 * D) + 769 * Sin(2 * L) - 668 *
    Sin(LS) - 412 * Sin(2 * F) - 212 * Sin(2 * L - 2 * D) - 206 * Sin(L + LS - 2 * D) + 192 *
    Sin(L + 2 * D) - 165 * Sin(LS - 2 * D) - 125 * Sin(D) - 110 * Sin(L + LS) + 148 * Sin(L - LS) -
    55 * Sin(2 * F - 2 * D);
  S := F + (DL + 412 * Sin(2 * F) + 541 * Sin(LS)) / ARC;
  H := F - 2 * D;
  N := -526 * Sin(H) + 44 * Sin(L + H) - 31 * Sin(-L + H) - 23 * Sin(LS + H) + 11 * Sin(-LS + H) -
    25 * Sin(-2 * L + F) + 21 * Sin(-L + F);
  L_MOON := P2 * Frac(L0 + DL / 1296E3); (* in rad *)
  B_MOON := (18520.0 * Sin(S) + N) / ARC; (* in rad *)
  (* equatorial coordinates *)
  CB := Cos(B_MOON);
  X := CB * Cos(L_MOON);
  V := CB * Sin(L_MOON);
  W := Sin(B_MOON);
  Y := COSEPS * V - SINEPS * W;
  Z := SINEPS * V + COSEPS * W;
  RHO := Sqrt(1.0 - Z * Z);
  Dec := (360.0 / P2) * ArcTan(Z / RHO);
  Ra := (48.0 / P2) * ArcTan(Y / (X + RHO));
  if Ra < 0 then
    Ra := Ra + 24.0;
end;

(*-----------------------------------------------------------------------*)

procedure Moon(T: Double; var LAMBDA, BETA, R: Double);

const
  PI2 = 6.283185308; (* 2*pi;  pi=3.141592654... *)
  ARC = 206264.81; (* 3600*180/pi = arcsec per radian *)

var
  DGAM, FAC: Double;
  DLAM, N, GAM1C, SINPI: Double;
  L0, L, LS, F, D, S: Double;
  DL0, DL, DLS, DF, DD, DS: Double;
  CO, SI: array [-6 .. 6, 1 .. 4] of Double;

  // fractional part of a number; with several compilers it may be
  // necessary to replace TRUNC by LongTrunc or INT if T<-24!
  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1;
    Frac := X
  end;

// calculate c=cos(a1+a2) and s=sin(a1+a2) from the addition theorems
// for c1=cos(a1), s1=sin(a1), c2=cos(a2) and s2=sin(a2)
  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

// calculate sin(phi); phi in units of 1 revolution = 360 degrees
  (*sub*)function Sine(PHI: Double): Double;
  begin
    Sine := Sin(PI2 * Frac(PHI));
  end;

// calculate long-periodic changes of the mean elements
// l,l',F,D and L0 as well as dgamma
  (*sub*)procedure Long_Periodic(T: Double; var DL0, DL, DLS, DF, DD, DGAM: Double);
  var
    S1, S2, S3, S4, S5, S6, S7: Double;
  begin
    S1 := Sine(0.19833 + 0.05611 * T);
    S2 := Sine(0.27869 + 0.04508 * T);
    S3 := Sine(0.16827 - 0.36903 * T);
    S4 := Sine(0.34734 - 5.37261 * T);
    S5 := Sine(0.10498 - 5.37899 * T);
    S6 := Sine(0.42681 - 0.41855 * T);
    S7 := Sine(0.14943 - 5.37511 * T);
    DL0 := 0.84 * S1 + 0.31 * S2 + 14.27 * S3 + 7.26 * S4 + 0.28 * S5 + 0.24 * S6;
    DL := 2.94 * S1 + 0.31 * S2 + 14.27 * S3 + 9.34 * S4 + 1.12 * S5 + 0.83 * S6;
    DLS := -6.40 * S1 - 1.89 * S6;
    DF := 0.21 * S1 + 0.31 * S2 + 14.27 * S3 - 88.70 * S4 - 15.30 * S5 + 0.24 * S6 - 1.86 * S7;
    DD := DL0 - DLS;
    DGAM := -3332E-9 * Sine(0.59734 - 5.37261 * T) - 539E-9 * Sine(0.35498 - 5.37899 * T) - 64E-9 *
      Sine(0.39943 - 5.37511 * T);
  end;

(* Init: calculates the mean elements and their sine and cosine *)
(* l mean anomaly of the Moon     l' mean anomaly of the Sun *)
(* F mean distance from the node  D  mean elongation from the Sun *)

  (*sub*)procedure Init;
  var
    I, J, MAX: integer;
    T2, ARG, FAC: Double;
  begin
    MAX := 1;
    T2 := T * T;
    DLAM := 0;
    DS := 0;
    GAM1C := 0;
    SINPI := 3422.7000;
    Long_Periodic(T, DL0, DL, DLS, DF, DD, DGAM);
    L0 := PI2 * Frac(0.60643382 + 1336.85522467 * T - 0.00000313 * T2) + DL0 / ARC;
    L := PI2 * Frac(0.37489701 + 1325.55240982 * T + 0.00002565 * T2) + DL / ARC;
    LS := PI2 * Frac(0.99312619 + 99.99735956 * T - 0.00000044 * T2) + DLS / ARC;
    F := PI2 * Frac(0.25909118 + 1342.22782980 * T - 0.00000892 * T2) + DF / ARC;
    D := PI2 * Frac(0.82736186 + 1236.85308708 * T - 0.00000397 * T2) + DD / ARC;
    for I := 1 to 4 do
    begin
      case I of
        1:
          begin
            ARG := L;
            MAX := 4;
            FAC := 1.000002208;
          end;
        2:
          begin
            ARG := LS;
            MAX := 3;
            FAC := 0.997504612 - 0.002495388 * T;
          end;
        3:
          begin
            ARG := F;
            MAX := 4;
            FAC := 1.000002708 + 139.978 * DGAM;
          end;
        4:
          begin
            ARG := D;
            MAX := 6;
            FAC := 1.0;
          end;
      end;
      CO[0, I] := 1.0;
      CO[1, I] := Cos(ARG) * FAC;
      SI[0, I] := 0.0;
      SI[1, I] := Sin(ARG) * FAC;
      for J := 2 to MAX do
        AddThe(CO[J - 1, I], SI[J - 1, I], CO[1, I], SI[1, I], CO[J, I], SI[J, I]);
      for J := 1 to MAX do
      begin
        CO[-J, I] := CO[J, I];
        SI[-J, I] := -SI[J, I];
      end;
    end;
  end;

// TERM calculates
// X=cos(p*arg1+q*arg2+r*arg3+s*arg4) and
// Y=sin(p*arg1+q*arg2+r*arg3+s*arg4)
  (*sub*)procedure Term(P, Q, R, S: integer; var X, Y: Double);
  var
    I: array [1 .. 4] of integer;
    K: integer;
  begin
    I[1] := P;
    I[2] := Q;
    I[3] := R;
    I[4] := S;
    X := 1.0;
    Y := 0.0;
    for K := 1 to 4 do
      if (I[K] <> 0) then
        AddThe(X, Y, CO[I[K], K], SI[I[K], K], X, Y);
  end;

  (*sub*)procedure AddSol(COEFFL, COEFFS, COEFFG, COEFFP: Double; P, Q, R, S: integer);
  var
    X, Y: Double;
  begin
    Term(P, Q, R, S, X, Y);
    DLAM := DLAM + COEFFL * Y;
    DS := DS + COEFFS * Y;
    GAM1C := GAM1C + COEFFG * X;
    SINPI := SINPI + COEFFP * X;
  end;

  (*sub*)procedure Solar1;
  begin
    AddSol(13.902, 14.06, -0.001, 0.2607, 0, 0, 0, 4);
    AddSol(0.403, -4.01, +0.394, 0.0023, 0, 0, 0, 3);
    AddSol(2369.912, 2373.36, +0.601, 28.2333, 0, 0, 0, 2);
    AddSol(-125.154, -112.79, -0.725, -0.9781, 0, 0, 0, 1);
    AddSol(1.979, 6.98, -0.445, 0.0433, 1, 0, 0, 4);
    AddSol(191.953, 192.72, +0.029, 3.0861, 1, 0, 0, 2);
    AddSol(-8.466, -13.51, +0.455, -0.1093, 1, 0, 0, 1);
    AddSol(22639.500, 22609.07, +0.079, 186.5398, 1, 0, 0, 0);
    AddSol(18.609, 3.59, -0.094, 0.0118, 1, 0, 0, -1);
    AddSol(-4586.465, -4578.13, -0.077, 34.3117, 1, 0, 0, -2);
    AddSol(+3.215, 5.44, +0.192, -0.0386, 1, 0, 0, -3);
    AddSol(-38.428, -38.64, +0.001, 0.6008, 1, 0, 0, -4);
    AddSol(-0.393, -1.43, -0.092, 0.0086, 1, 0, 0, -6);
    AddSol(-0.289, -1.59, +0.123, -0.0053, 0, 1, 0, 4);
    AddSol(-24.420, -25.10, +0.040, -0.3000, 0, 1, 0, 2);
    AddSol(18.023, 17.93, +0.007, 0.1494, 0, 1, 0, 1);
    AddSol(-668.146, -126.98, -1.302, -0.3997, 0, 1, 0, 0);
    AddSol(0.560, 0.32, -0.001, -0.0037, 0, 1, 0, -1);
    AddSol(-165.145, -165.06, +0.054, 1.9178, 0, 1, 0, -2);
    AddSol(-1.877, -6.46, -0.416, 0.0339, 0, 1, 0, -4);
    AddSol(0.213, 1.02, -0.074, 0.0054, 2, 0, 0, 4);
    AddSol(14.387, 14.78, -0.017, 0.2833, 2, 0, 0, 2);
    AddSol(-0.586, -1.20, +0.054, -0.0100, 2, 0, 0, 1);
    AddSol(769.016, 767.96, +0.107, 10.1657, 2, 0, 0, 0);
    AddSol(+1.750, 2.01, -0.018, 0.0155, 2, 0, 0, -1);
    AddSol(-211.656, -152.53, +5.679, -0.3039, 2, 0, 0, -2);
    AddSol(+1.225, 0.91, -0.030, -0.0088, 2, 0, 0, -3);
    AddSol(-30.773, -34.07, -0.308, 0.3722, 2, 0, 0, -4);
    AddSol(-0.570, -1.40, -0.074, 0.0109, 2, 0, 0, -6);
    AddSol(-2.921, -11.75, +0.787, -0.0484, 1, 1, 0, 2);
    AddSol(+1.267, 1.52, -0.022, 0.0164, 1, 1, 0, 1);
    AddSol(-109.673, -115.18, +0.461, -0.9490, 1, 1, 0, 0);
    AddSol(-205.962, -182.36, +2.056, +1.4437, 1, 1, 0, -2);
    AddSol(0.233, 0.36, 0.012, -0.0025, 1, 1, 0, -3);
    AddSol(-4.391, -9.66, -0.471, 0.0673, 1, 1, 0, -4);
  end;

  (*sub*)procedure Solar2;
  begin
    AddSol(0.283, 1.53, -0.111, +0.0060, 1, -1, 0, +4);
    AddSol(14.577, 31.70, -1.540, +0.2302, 1, -1, 0, 2);
    AddSol(147.687, 138.76, +0.679, +1.1528, 1, -1, 0, 0);
    AddSol(-1.089, 0.55, +0.021, 0.0, 1, -1, 0, -1);
    AddSol(28.475, 23.59, -0.443, -0.2257, 1, -1, 0, -2);
    AddSol(-0.276, -0.38, -0.006, -0.0036, 1, -1, 0, -3);
    AddSol(0.636, 2.27, +0.146, -0.0102, 1, -1, 0, -4);
    AddSol(-0.189, -1.68, +0.131, -0.0028, 0, 2, 0, 2);
    AddSol(-7.486, -0.66, -0.037, -0.0086, 0, 2, 0, 0);
    AddSol(-8.096, -16.35, -0.740, 0.0918, 0, 2, 0, -2);
    AddSol(-5.741, -0.04, 0.0, -0.0009, 0, 0, 2, 2);
    AddSol(0.255, 0.0, 0.0, 0.0, 0, 0, 2, 1);
    AddSol(-411.608, -0.20, 0.0, -0.0124, 0, 0, 2, 0);
    AddSol(0.584, 0.84, 0.0, +0.0071, 0, 0, 2, -1);
    AddSol(-55.173, -52.14, 0.0, -0.1052, 0, 0, 2, -2);
    AddSol(0.254, 0.25, 0.0, -0.0017, 0, 0, 2, -3);
    AddSol(+0.025, -1.67, 0.0, +0.0031, 0, 0, 2, -4);
    AddSol(1.060, 2.96, -0.166, 0.0243, 3, 0, 0, +2);
    AddSol(36.124, 50.64, -1.300, 0.6215, 3, 0, 0, 0);
    AddSol(-13.193, -16.40, +0.258, -0.1187, 3, 0, 0, -2);
    AddSol(-1.187, -0.74, +0.042, 0.0074, 3, 0, 0, -4);
    AddSol(-0.293, -0.31, -0.002, 0.0046, 3, 0, 0, -6);
    AddSol(-0.290, -1.45, +0.116, -0.0051, 2, 1, 0, 2);
    AddSol(-7.649, -10.56, +0.259, -0.1038, 2, 1, 0, 0);
    AddSol(-8.627, -7.59, +0.078, -0.0192, 2, 1, 0, -2);
    AddSol(-2.740, -2.54, +0.022, 0.0324, 2, 1, 0, -4);
    AddSol(1.181, 3.32, -0.212, 0.0213, 2, -1, 0, +2);
    AddSol(9.703, 11.67, -0.151, 0.1268, 2, -1, 0, 0);
    AddSol(-0.352, -0.37, +0.001, -0.0028, 2, -1, 0, -1);
    AddSol(-2.494, -1.17, -0.003, -0.0017, 2, -1, 0, -2);
    AddSol(0.360, 0.20, -0.012, -0.0043, 2, -1, 0, -4);
    AddSol(-1.167, -1.25, +0.008, -0.0106, 1, 2, 0, 0);
    AddSol(-7.412, -6.12, +0.117, 0.0484, 1, 2, 0, -2);
    AddSol(-0.311, -0.65, -0.032, 0.0044, 1, 2, 0, -4);
    AddSol(+0.757, 1.82, -0.105, 0.0112, 1, -2, 0, 2);
    AddSol(+2.580, 2.32, +0.027, 0.0196, 1, -2, 0, 0);
    AddSol(+2.533, 2.40, -0.014, -0.0212, 1, -2, 0, -2);
    AddSol(-0.344, -0.57, -0.025, +0.0036, 0, 3, 0, -2);
    AddSol(-0.992, -0.02, 0.0, 0.0, 1, 0, 2, 2);
    AddSol(-45.099, -0.02, 0.0, -0.0010, 1, 0, 2, 0);
    AddSol(-0.179, -9.52, 0.0, -0.0833, 1, 0, 2, -2);
    AddSol(-0.301, -0.33, 0.0, 0.0014, 1, 0, 2, -4);
    AddSol(-6.382, -3.37, 0.0, -0.0481, 1, 0, -2, 2);
    AddSol(39.528, 85.13, 0.0, -0.7136, 1, 0, -2, 0);
    AddSol(9.366, 0.71, 0.0, -0.0112, 1, 0, -2, -2);
    AddSol(0.202, 0.02, 0.0, 0.0, 1, 0, -2, -4);
  end;

  (*sub*)procedure Solar3;
  begin
    AddSol(0.415, 0.10, 0.0, 0.0013, 0, 1, 2, 0);
    AddSol(-2.152, -2.26, 0.0, -0.0066, 0, 1, 2, -2);
    AddSol(-1.440, -1.30, 0.0, +0.0014, 0, 1, -2, 2);
    AddSol(0.384, -0.04, 0.0, 0.0, 0, 1, -2, -2);
    AddSol(+1.938, +3.60, -0.145, +0.0401, 4, 0, 0, 0);
    AddSol(-0.952, -1.58, +0.052, -0.0130, 4, 0, 0, -2);
    AddSol(-0.551, -0.94, +0.032, -0.0097, 3, 1, 0, 0);
    AddSol(-0.482, -0.57, +0.005, -0.0045, 3, 1, 0, -2);
    AddSol(0.681, 0.96, -0.026, 0.0115, 3, -1, 0, 0);
    AddSol(-0.297, -0.27, 0.002, -0.0009, 2, 2, 0, -2);
    AddSol(0.254, +0.21, -0.003, 0.0, 2, -2, 0, -2);
    AddSol(-0.250, -0.22, 0.004, 0.0014, 1, 3, 0, -2);
    AddSol(-3.996, 0.0, 0.0, +0.0004, 2, 0, 2, 0);
    AddSol(0.557, -0.75, 0.0, -0.0090, 2, 0, 2, -2);
    AddSol(-0.459, -0.38, 0.0, -0.0053, 2, 0, -2, 2);
    AddSol(-1.298, 0.74, 0.0, +0.0004, 2, 0, -2, 0);
    AddSol(0.538, 1.14, 0.0, -0.0141, 2, 0, -2, -2);
    AddSol(0.263, 0.02, 0.0, 0.0, 1, 1, 2, 0);
    AddSol(0.426, +0.07, 0.0, -0.0006, 1, 1, -2, -2);
    AddSol(-0.304, +0.03, 0.0, +0.0003, 1, -1, 2, 0);
    AddSol(-0.372, -0.19, 0.0, -0.0027, 1, -1, -2, 2);
    AddSol(+0.418, 0.0, 0.0, 0.0, 0, 0, 4, 0);
    AddSol(-0.330, -0.04, 0.0, 0.0, 3, 0, 2, 0);
  end;

// part N of the perturbations of ecliptic latitude
  (*sub*)procedure SolarN(var N: Double);
  var
    X, Y: Double;
    (*sub*)(*sub*)procedure AddN(COEFFN: Double; P, Q, R, S: integer);
    begin
      Term(P, Q, R, S, X, Y);
      N := N + COEFFN * Y
    end;

  begin
    N := 0.0;
    AddN(-526.069, 0, 0, 1, -2);
    AddN(-3.352, 0, 0, 1, -4);
    AddN(+44.297, +1, 0, 1, -2);
    AddN(-6.000, +1, 0, 1, -4);
    AddN(+20.599, -1, 0, 1, 0);
    AddN(-30.598, -1, 0, 1, -2);
    AddN(-24.649, -2, 0, 1, 0);
    AddN(-2.000, -2, 0, 1, -2);
    AddN(-22.571, 0, +1, 1, -2);
    AddN(+10.985, 0, -1, 1, -2);
  end;

(* perturbations of ecliptic latitude by Venus and Jupiter *)
  (*sub*)procedure Planetary(var DLAM: Double);
  begin
    DLAM := DLAM + 0.82 * Sine(0.7736 - 62.5512 * T) + 0.31 * Sine(0.0466 - 125.1025 * T) + 0.35 *
      Sine(0.5785 - 25.1042 * T) + 0.66 * Sine(0.4591 + 1335.8075 * T) + 0.64 *
      Sine(0.3130 - 91.5680 * T) + 1.14 * Sine(0.1480 + 1331.2898 * T) + 0.21 *
      Sine(0.5918 + 1056.5859 * T) + 0.44 * Sine(0.5784 + 1322.8595 * T) + 0.24 *
      Sine(0.2275 - 5.7374 * T) + 0.28 * Sine(0.2965 + 2.6929 * T) + 0.33 *
      Sine(0.3132 + 6.3368 * T);
  end;

begin

  Init;
  Solar1;
  Solar2;
  Solar3;
  SolarN(N);
  Planetary(DLAM);

  LAMBDA := 360.0 * Frac((L0 + DLAM / ARC) / PI2);

  S := F + DS / ARC;
  FAC := 1.000002708 + 139.978 * DGAM;
  BETA := (FAC * (18518.511 + 1.189 + GAM1C) * Sin(S) - 6.24 * Sin(3 * S) + N) / 3600.0;

  SINPI := SINPI * 0.999953253;
  R := ARC / SINPI;
end;

(* ----------------------------------------------------------------------- *)
procedure MoonEqu(T: Double; var Ra, Dec, R: Double);
var
  L, B, X, Y, Z: Double;
begin
  Moon(T, L, B, R); // ecliptic coordinates (mean equinox
  Cart(R, B, L, X, Y, Z); // of date)
  Ecl2Equ(T, X, Y, Z); // transform into equatorial coordinates
  NutEqu(T, X, Y, Z); // nutation
  Polar(X, Y, Z, R, Dec, Ra);
end;

(* ???
procedure T_Fit_Moon(TA, TB: Double; N: integer; var RA_POLY, DE_POLY, R_POLY: TPolynomCheb);
begin
  T_Fit_LBR({MoonEqu,?} TA, TB, N, RA_POLY, DE_POLY, R_POLY);
end;
*)

(*-----------------------------------------------------------------------*)
procedure T_Fit_Moon(TA, TB: Double; N: integer; var RA_POLY, DE_POLY, R_POLY: TPolynomCheb);
const
  NDIM = 27;
var
  I, J, K: integer;
  FAC, BPA, BMA, PHI: Double;
  T, H, Ra, DE, R: array [0 .. NDIM] of Double;
begin
  if (NDIM < 2 * MAX_TP_DEG + 1) then
    writeln(' NDIM too small in T_Fit_Moon');
  if (N > MAX_TP_DEG) then
    writeln(' N too large in T_Fit_Moon');
  RA_POLY.M := N;
  DE_POLY.M := N;
  R_POLY.M := N;
  RA_POLY.A := TA;
  DE_POLY.A := TA;
  R_POLY.A := TA;
  RA_POLY.B := TB;
  DE_POLY.B := TB;
  R_POLY.B := TB;
  BMA := (TB - TA) / 2.0;
  BPA := (TB + TA) / 2.0;
  FAC := 2.0 / (N + 1);
  PHI := Pi / (2 * N + 2); (* h(k)=cos(pi*k/N/2) *)
  H[0] := 1.0;
  H[1] := Cos(PHI);
  for I := 2 to (2 * N + 1) do
    H[I] := 2 * H[1] * H[I - 1] - H[I - 2];
  for K := 1 to N + 1 do
    T[K] := H[2 * K - 1] * BMA + BPA; (* subdivision points *)
  for K := 1 to N + 1 do
    MoonEqu(T[K], Ra[K], DE[K], R[K]);
  for K := 2 to N + 1 do (* make Ra continuous *)
    if (Ra[K - 1] < Ra[K]) then
      Ra[K] := Ra[K] - 360.0; (* in [-360,+360] deg *)
  for J := 0 to N do (* calculate Chebyshev *)
  begin (* coefficients C[J] *)
    PHI := Pi * J / (2 * N + 2);
    H[1] := Cos(PHI);
    for I := 2 to (2 * N + 1) do
      H[I] := 2 * H[1] * H[I - 1] - H[I - 2];
    RA_POLY.C[J] := 0.0;
    DE_POLY.C[J] := 0.0;
    R_POLY.C[J] := 0.0;
    for K := 1 to N + 1 do
    begin
      RA_POLY.C[J] := RA_POLY.C[J] + H[2 * K - 1] * Ra[K];
      DE_POLY.C[J] := DE_POLY.C[J] + H[2 * K - 1] * DE[K];
      R_POLY.C[J] := R_POLY.C[J] + H[2 * K - 1] * R[K];
    end;
    RA_POLY.C[J] := RA_POLY.C[J] * FAC;
    DE_POLY.C[J] := DE_POLY.C[J] * FAC;
    R_POLY.C[J] := R_POLY.C[J] * FAC;
  end;
end;



end.
