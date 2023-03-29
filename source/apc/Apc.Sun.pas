unit Apc.Sun;

interface

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric;

(*-----------------------------------------------------------------------*)
(* SunPos: ecliptic coordinates L,B,R (in deg and AU) of the             *)
(*         Sun referred to the mean equinox of date                      *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure SunPos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* SunEqu: apparent equatorial coordinates of the Sun                    *)
(*         (right ascension Ra, declination Dec in deg, R in AU)         *)
(*         T in Julian centuries since J2000 (T=(JD-2451545.0)/36525)    *)
(*-----------------------------------------------------------------------*)
procedure SunEqu(T: Double; var Ra, Dec, R: Double);

(*-----------------------------------------------------------------------*)
(* MiniSun: low precision solar coordinates (approx. 1')                *)
(*           T  : time in Julian centuries since J2000                   *)
(*                ( T=(JD-2451545)/36525 )                               *)
(*           Ra : right ascension (in h; equinox of date)                *)
(*           Dec: declination (in deg; equinox of date)                  *)
(*-----------------------------------------------------------------------*)
procedure MiniSun(T: Double; var Ra, Dec: Double);

(*-----------------------------------------------------------------------*)
(* T_Fit_Sun: approximates the equatorial coordinates of the             *)
(*            Sun by Chebyshev expansions for a given period of time     *)
(*                                                                       *)
(*  TA     : first date (in Julian centuries since J2000)                *)
(*  TB     : last date ( TB < TA + 1 year )                              *)
(*  N      : highest order                                               *)
(*  RA_POLY: coefficients for right ascension                            *)
(*  DE_POLY: coefficients for declination                                *)
(*  R_POLY : coefficients for geocentric distance                        *)
(*-----------------------------------------------------------------------*)
procedure T_Fit_Sun(TA, TB: Double; N: integer; var RA_POLY, DE_POLY, R_POLY: TPolynomCheb);

//=======================================================================
implementation
//=======================================================================

(*-----------------------------------------------------------------------*)

procedure MiniSun(T: Double; var Ra, Dec: Double);
const
  P2 = 6.283185307;
  COSEPS = 0.91748;
  SINEPS = 0.39778;
var
  L, M, DL, SL, X, Y, Z, RHO: Double;
  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1;
    Frac := X
  end;

begin
  M := P2 * Frac(0.993133 + 99.997361 * T);
  DL := 6893.0 * Sin(M) + 72.0 * Sin(2 * M);
  L := P2 * Frac(0.7859453 + M / P2 + (6191.2 * T + DL) / 1296E3);
  SL := Sin(L);
  X := Cos(L);
  Y := COSEPS * SL;
  Z := SINEPS * SL;
  RHO := Sqrt(1.0 - Z * Z);
  Dec := (360.0 / P2) * ArcTan(Z / RHO);
  Ra := (48.0 / P2) * ArcTan(Y / (X + RHO));
  if (Ra < 0) then
    Ra := Ra + 24.0;
end;  // MiniSun

(* ----------------------------------------------------------------------- *)

procedure SunPos(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C3, S3: array [-1 .. 7] of Double;
  C, S: array [-8 .. 0] of Double;
  M2, M3, M4, M5, M6: Double;
  D, A, UU: Double;
  U, V, DL, DR, DB: Double;
  I: integer;

  (*sub*)function Frac(X: Double): Double;
  (* for some compilers Trunc has to be replaced by LongTrunc *)
  (* or INT if the routine is used with T<-24 *)
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I1, I, IT: integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C3[I1], S3[I1], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertVenus; (* Keplerian terms and perturbations by Venus *)
  var
    I: integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M2);
    S[-1] := -Sin(M2);
    for I := -1 downto -5 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, 0, 0, -0.22, 6892.76, -16707.37, -0.54, 0.00, 0.00);
    Term(1, 0, 1, -0.06, -17.35, 42.04, -0.15, 0.00, 0.00);
    Term(1, 0, 2, -0.01, -0.05, 0.13, -0.02, 0.00, 0.00);
    Term(2, 0, 0, 0.00, 71.98, -139.57, 0.00, 0.00, 0.00);
    Term(2, 0, 1, 0.00, -0.36, 0.70, 0.00, 0.00, 0.00);
    Term(3, 0, 0, 0.00, 1.04, -1.75, 0.00, 0.00, 0.00);
    Term(0, -1, 0, 0.03, -0.07, -0.16, -0.07, 0.02, -0.02);
    Term(1, -1, 0, 2.35, -4.23, -4.75, -2.64, 0.00, 0.00);
    Term(1, -2, 0, -0.10, 0.06, 0.12, 0.20, 0.02, 0.00);
    Term(2, -1, 0, -0.06, -0.03, 0.20, -0.01, 0.01, -0.09);
    Term(2, -2, 0, -4.70, 2.90, 8.28, 13.42, 0.01, -0.01);
    Term(3, -2, 0, 1.80, -1.74, -1.44, -1.57, 0.04, -0.06);
    Term(3, -3, 0, -0.67, 0.03, 0.11, 2.43, 0.01, 0.00);
    Term(4, -2, 0, 0.03, -0.03, 0.10, 0.09, 0.01, -0.01);
    Term(4, -3, 0, 1.51, -0.40, -0.88, -3.36, 0.18, -0.10);
    Term(4, -4, 0, -0.19, -0.09, -0.38, 0.77, 0.00, 0.00);
    Term(5, -3, 0, 0.76, -0.68, 0.30, 0.37, 0.01, 0.00);
    Term(5, -4, 0, -0.14, -0.04, -0.11, 0.43, -0.03, 0.00);
    Term(5, -5, 0, -0.05, -0.07, -0.31, 0.21, 0.00, 0.00);
    Term(6, -4, 0, 0.15, -0.04, -0.06, -0.21, 0.01, 0.00);
    Term(6, -5, 0, -0.03, -0.03, -0.09, 0.09, -0.01, 0.00);
    Term(6, -6, 0, 0.00, -0.04, -0.18, 0.02, 0.00, 0.00);
    Term(7, -5, 0, -0.12, -0.03, -0.08, 0.31, -0.02, -0.01);
  end;

  (*sub*)procedure PertMars; (* perturbations by Mars *)
  var
    I: integer;
  begin
    C[-1] := Cos(M4);
    S[-1] := -Sin(M4);
    for I := -1 downto -7 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, -1, 0, -0.22, 0.17, -0.21, -0.27, 0.00, 0.00);
    Term(1, -2, 0, -1.66, 0.62, 0.16, 0.28, 0.00, 0.00);
    Term(2, -2, 0, 1.96, 0.57, -1.32, 4.55, 0.00, 0.01);
    Term(2, -3, 0, 0.40, 0.15, -0.17, 0.46, 0.00, 0.00);
    Term(2, -4, 0, 0.53, 0.26, 0.09, -0.22, 0.00, 0.00);
    Term(3, -3, 0, 0.05, 0.12, -0.35, 0.15, 0.00, 0.00);
    Term(3, -4, 0, -0.13, -0.48, 1.06, -0.29, 0.01, 0.00);
    Term(3, -5, 0, -0.04, -0.20, 0.20, -0.04, 0.00, 0.00);
    Term(4, -4, 0, 0.00, -0.03, 0.10, 0.04, 0.00, 0.00);
    Term(4, -5, 0, 0.05, -0.07, 0.20, 0.14, 0.00, 0.00);
    Term(4, -6, 0, -0.10, 0.11, -0.23, -0.22, 0.00, 0.00);
    Term(5, -7, 0, -0.05, 0.00, 0.01, -0.14, 0.00, 0.00);
    Term(5, -8, 0, 0.05, 0.01, -0.02, 0.10, 0.00, 0.00);
  end;

  (*sub*)procedure PertJupiter; (* perturbations by Jupiter *)
  var
    I: integer;
  begin
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    for I := -1 downto -3 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    TERM(-1, -1, 0, 0.01, 0.07, 0.18, -0.02, 0.00, -0.02);
    Term(0, -1, 0, -0.31, 2.58, 0.52, 0.34, 0.02, 0.00);
    Term(1, -1, 0, -7.21, -0.06, 0.13, -16.27, 0.00, -0.02);
    Term(1, -2, 0, -0.54, -1.52, 3.09, -1.12, 0.01, -0.17);
    Term(1, -3, 0, -0.03, -0.21, 0.38, -0.06, 0.00, -0.02);
    Term(2, -1, 0, -0.16, 0.05, -0.18, -0.31, 0.01, 0.00);
    Term(2, -2, 0, 0.14, -2.73, 9.23, 0.48, 0.00, 0.00);
    Term(2, -3, 0, 0.07, -0.55, 1.83, 0.25, 0.01, 0.00);
    Term(2, -4, 0, 0.02, -0.08, 0.25, 0.06, 0.00, 0.00);
    Term(3, -2, 0, 0.01, -0.07, 0.16, 0.04, 0.00, 0.00);
    Term(3, -3, 0, -0.16, -0.03, 0.08, -0.64, 0.00, 0.00);
    Term(3, -4, 0, -0.04, -0.01, 0.03, -0.17, 0.00, 0.00);
  end;

  (*sub*)procedure PertSaturn; (* perturbations by Saturn *)
  begin
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(0, -1, 0, 0.00, 0.32, 0.01, 0.00, 0.00, 0.00);
    Term(1, -1, 0, -0.08, -0.41, 0.97, -0.18, 0.00, -0.01);
    Term(1, -2, 0, 0.04, 0.10, -0.23, 0.10, 0.00, 0.00);
    Term(2, -2, 0, 0.04, 0.10, -0.35, 0.13, 0.00, 0.00);
  end;

  (*sub*)procedure PertMoon; (* difference between the Earth-Moon *)
  begin (* barycenter and the center of the Earth *)
    DL := DL + 6.45 * Sin(D) - 0.42 * Sin(D - A) + 0.18 * Sin(D + A) + 0.17 * Sin(D - M3) - 0.06 *
      Sin(D + M3);
    DR := DR + 30.76 * Cos(D) - 3.06 * Cos(D - A) + 0.85 * Cos(D + A) - 0.58 * Cos(D + M3) + 0.57 *
      Cos(D - M3);
    DB := DB + 0.576 * Sin(UU);
  end;

begin (* SunPos *)
  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M2 := P2 * Frac(0.1387306 + 162.5485917 * T);
  M3 := P2 * Frac(0.9931266 + 99.9973604 * T);
  M4 := P2 * Frac(0.0543250 + 53.1666028 * T);
  M5 := P2 * Frac(0.0551750 + 8.4293972 * T);
  M6 := P2 * Frac(0.8816500 + 3.3938722 * T);
  D := P2 * Frac(0.8274 + 1236.8531 * T);
  A := P2 * Frac(0.3749 + 1325.5524 * T);
  UU := P2 * Frac(0.2591 + 1342.2278 * T);
  C3[0] := 1.0;
  S3[0] := 0.0;
  C3[1] := Cos(M3);
  S3[1] := Sin(M3);
  C3[-1] := C3[1];
  S3[-1] := -S3[1];
  for I := 2 to 7 do
    AddThe(C3[I - 1], S3[I - 1], C3[1], S3[1], C3[I], S3[I]);
  PertVenus;
  PertMars;
  PertJupiter;
  PertSaturn;
  PertMoon;
  DL := DL + 6.40 * Sin(P2 * (0.6983 + 0.0561 * T)) + 1.87 * Sin(P2 * (0.5764 + 0.4174 * T)) + 0.27
    * Sin(P2 * (0.4189 + 0.3306 * T)) + 0.20 * Sin(P2 * (0.3581 + 2.4814 * T));
  L := 360.0 * Frac(0.7859453 + M3 / P2 + ((6191.2 + 1.1 * T) * T + DL) / 1296.0E3);
  R := 1.0001398 - 0.0000007 * T + DR * 1E-6;
  B := DB / 3600.0;

end; (* SunPos *)

(* ----------------------------------------------------------------------- *)

procedure SunEqu(T: Double; var Ra, Dec, R: Double);
var
  DT, L, B, X, Y, Z: Double;
begin
  DT := (8.32 / 1440.0) / 36525.0; (* light-time correction of 8.32 min *)
  SunPos(T - DT, L, B, R); (* geocentric ecliptic coordinates *)
  Cart(R, B, L, X, Y, Z); (* cartesian ecliptic coordinates *)
  Ecl2Equ(T, X, Y, Z); (* equatorial ecliptic coordinates *)
  NutEqu(T, X, Y, Z); (* correction for nutation *)
  Polar(X, Y, Z, R, Dec, Ra); (* spherical coordinates; *)
end; (* true equinox of date *)

(*  ???
procedure T_Fit_Sun(TA, TB: Double; N: integer; var RA_POLY, DE_POLY, R_POLY: TPolynomCheb);
begin
  T_Fit_LBR({SunEqu,?} TA, TB, N, RA_POLY, DE_POLY, R_POLY);
end;
*)

(*-----------------------------------------------------------------------------------*)

procedure T_Fit_Sun(TA, TB: Double; N: integer; var RA_POLY, DE_POLY, R_POLY: TPolynomCheb);
const
  Pi = 3.1415926535898;
  NDIM = 27;
var
  I, J, K: integer;
  FAC, BPA, BMA, PHI: Double;
  T, H, Ra, DE, R: array [0 .. NDIM] of Double;
begin
  if (NDIM < 2 * MAX_TP_DEG + 1) then
    writeln(' NDIM too small in T_Fit_Sun');
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
    SunEqu(T[K], Ra[K], DE[K], R[K]);
  for K := 2 to N + 1 do (* make Ra continuous *)
    if (Ra[K - 1] < Ra[K]) then
      Ra[K] := Ra[K] - 360.0; (* in [-360,+360] deg *)
  for J := 0 to N do (* calculate Chebyshev *)
  begin (* coefficients C(j) *)
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
end;  // T_Fit_Sun


end.
