unit Astro.Matlib;

interface

uses
  System.Math;

(* vectors and matrices *)

type
  INDEX = (X, Y, Z);
  VECTOR = array [INDEX] of Double;
  REAL3 = array [1 .. 3] of Double;
  Double33 = array [1 .. 3] of REAL3;
  MAT3X = array [1 .. 3] of VECTOR;

  (* Chebyshev polynomials *)

const
  MAX_TP_DEG = 13; (* maximum degree *)

type
  TPolynom = record (* Chebyshev polynomial *)
    M: integer; (* degree *)
    A, B: Double; (* interval *)
    C: array [0 .. MAX_TP_DEG] of Double; (* coefficients *)
  end;

  (* Vector and matrix for least squares systems: *)
  (* minimum dimensions for use with program FOTO are given below; *)
  (* they may be increased as needed by the user. *)

type
  LSQVEC = array [1 .. 3] of Double;
  LSQMAT = array [1 .. 30, 1 .. 5] of Double;


(* functions and procedures *)

(* -------------------------------------------------------------------- *)
(* ACS: arccosine function (degrees) *)
(* -------------------------------------------------------------------- *)
function ACS(X: double): double;
(* -------------------------------------------------------------------- *)
(* ASN: arcsine function (degrees) *)
(* -------------------------------------------------------------------- *)
function ASN(X: double): double;
(*-----------------------------------------------------------------------*)
(* ATN: arctangent function (degrees)                                    *)
(*-----------------------------------------------------------------------*)
function ATN(X: double): double;
(*-----------------------------------------------------------------------*)
(* ATN2: arctangent of y/x for two arguments                             *)
(*       (correct quadrant; -180 deg <= ATN2 <= +180 deg)                *)
(*-----------------------------------------------------------------------*)
function ATN2(Y,X: double): double;
(*-----------------------------------------------------------------------*)
(* CART: conversion of polar coordinates (r,theta,phi)                   *)
(*       into cartesian coordinates (x,y,z)                              *)
(*       (theta in [-90 deg,+90 deg]; phi in [-360 deg,+360 deg])        *)
(*-----------------------------------------------------------------------*)
procedure CART(R,THETA,PHI: Double; var X,Y,Z: Double);
(*-----------------------------------------------------------------------*)
(* CROSS: cross product of two vectors                                   *)
(*-----------------------------------------------------------------------*)
procedure CROSS(A,B:VECTOR;var C:VECTOR);
(*-----------------------------------------------------------------------*)
(* CS: cosine function (degrees)                                         *)
(*-----------------------------------------------------------------------*)
function CS(X: double): double;
(*-----------------------------------------------------------------------*)
(* CUBR: cube root                                                       *)
(*-----------------------------------------------------------------------*)
function CUBR(X: double): double;
(*-----------------------------------------------------------------------*)
(* DDD: conversion of degrees, minutes and seconds into                  *)
(*      degrees and fractions of a degree                                *)
(*-----------------------------------------------------------------------*)
procedure DDD(D, M: integer; S: Double; var DD: Double);
(*-----------------------------------------------------------------------*)
(* DMS: conversion of degrees and fractions of a degree                  *)
(*      into degrees, minutes and seconds                                *)
(*-----------------------------------------------------------------------*)
procedure DMS(DD: Double; var D, M: integer; var S: Double);
(*-----------------------------------------------------------------------*)
(* DOT: dot product of two vectors                                       *)
(*-----------------------------------------------------------------------*)
function DOT(A, B: VECTOR): Double;
(*-----------------------------------------------------------------------*)
(* LSQFIT:                                                               *)
(*   solution of an overdetermined system of linear equations            *)
(*   A[i,1]*s[1]+...A[i,m]*s[m] - A[i,m+1] = 0   (i=1,..,n)              *)
(*   according to the method of least squares using Givens rotations     *)
(*   A: matrix of coefficients                                           *)
(*   N: number of equations  (rows of A)                                 *)
(*   M: number of unknowns   (M+1=columns of A, M=elements of S)         *)
(*   S: solution vector                                                  *)
(*-----------------------------------------------------------------------*)
procedure LSQFIT(A: LSQMAT; N, M: integer; var S: LSQVEC);
(*-----------------------------------------------------------------------*)
(* NORM: magnitude of a vector                                           *)
(*-----------------------------------------------------------------------*)
function NORM(A: VECTOR): Double;
(*-----------------------------------------------------------------------*)
(* POLAR: conversion of cartesian coordinates (x,y,z)                    *)
(*        into polar coordinates (r,theta,phi)                           *)
(*        (theta in [-90 deg,+90 deg]; phi in [0 deg,+360 deg])          *)
(*-----------------------------------------------------------------------*)
procedure POLAR(X, Y, Z: Double; var R, THETA, PHI: Double);
(*-----------------------------------------------------------------------*)
(* QUAD: finds a parabola through 3 points                               *)
(*       (-1,Y_MINUS), (0,Y_0) und (1,Y_PLUS),                           *)
(*       that do not lie on a straight line.                             *)
(*                                                                       *)
(*      Y_MINUS,Y_0,Y_PLUS: three y-values                               *)
(*      XE,YE   : x and y of the extreme value of the parabola           *)
(*      ZERO1   : first root within [-1,+1] (for NZ=1,2)                 *)
(*      ZERO2   : second root within [-1,+1] (only for NZ=2)             *)
(*      NZ      : number of roots within the interval [-1,+1]            *)
(*-----------------------------------------------------------------------*)
procedure QUAD(Y_MINUS, Y_0, Y_PLUS: Double; var XE, YE, ZERO1, ZERO2: Double; var NZ: integer);
(*-----------------------------------------------------------------------*)
(* T_EVAL: evaluates the approximation of a function by Chebyshev        *)
(*         polynomials of maximum order F.M over the interval [F.A,F.B]  *)
(*  F : record containing the Chebyshev coefficients                     *)
(*  X : argument                                                         *)
(*-----------------------------------------------------------------------*)
 function T_EVAL(F: TPolynom; X: Double): Double;
(*-----------------------------------------------------------------------*)
(* T_FIT_LBR: expands lunar or planetary coordinates into series of      *)
(*            Chebyshev polynomials for longitude, latitude and radius   *)
(*            that are valid for a specified period of time              *)
(*                                                                       *)
(*  procedure POSITION: routine for calculating the coordinates L,B,R    *)
(*  TA      : first date of desired period of time                       *)
(*  TB      : last date                                                  *)
(*  N       : highest order of Chebyshev polynomials (N<=MAX_TP_DEG)     *)
(*  L_POLY  : coefficients for longitude                                 *)
(*  B_POLY  : coefficients for latitude                                  *)
(*  R_POLY  : coefficients for radius                                    *)
(*                                                                       *)
(* note:                                                                 *)
(*  . the interval [TA,TB] must be shorter than one revolution!          *)
(*  . the routine will only work for heliocentric planetary or           *)
(*    geocentric lunar but not for geocentric planetary coordinates!     *)
(*-----------------------------------------------------------------------*)
procedure T_FIT_LBR(TA, TB: Double; N: integer; var L_POLY, B_POLY, R_POLY: TPolynom);

(*-----------------------------------------------------------------------*)
(* SN: sine function (degrees)                                           *)
(*-----------------------------------------------------------------------*)
function SN(X: double): double;
(*-----------------------------------------------------------------------*)
(* TN: tangent function (degrees)                                        *)
(*-----------------------------------------------------------------------*)
function TN(X: Double): Double;


implementation

(* ----------------------------------------------------------------------- *)

function ACS(X: double): double;
const
  RAD = 0.0174532925199433;
  EPS = 1E-7;
  C = 90.0;
begin
  if abs(X) = 1.0 then
    ACS := C - X * C
  else if (abs(X) > EPS) then
    ACS := C - ArcTan(X / Sqrt((1.0 - X) * (1.0 + X))) / RAD
  else
    ACS := C - X / RAD;
end;

(* ----------------------------------------------------------------------- *)

function ASN(X: double): double;
const
  RAD = 0.0174532925199433;
  EPS = 1E-7;
begin
  if abs(X) = 1.0 then
    ASN := 90.0 * X
  else if (abs(X) > EPS) then
    ASN := ArcTan(X / Sqrt((1.0 - X) * (1.0 + X))) / RAD
  else
    ASN := X / RAD;
end;

(*-----------------------------------------------------------------------*)
function ATN(X: double): double;
const
  RAD = 0.0174532925199433;
begin
  ATN := ArcTan(X) / RAD
end;

(*-----------------------------------------------------------------------*)

function ATN2(Y, X: double): double;
const
  RAD = 0.0174532925199433;
var
  AX, AY, PHI: Double;
begin
  if (X = 0.0) and (Y = 0.0) then
    ATN2 := 0.0
  else
  begin
    AX := abs(X);
    AY := abs(Y);
    if (AX > AY) then
      PHI := ArcTan(AY / AX) / RAD
    else
      PHI := 90.0 - ArcTan(AX / AY) / RAD;
    if (X < 0.0) then
      PHI := 180.0 - PHI;
    if (Y < 0.0) then
      PHI := -PHI;
    ATN2 := PHI;
  end;
end;

(*-----------------------------------------------------------------------*)
procedure CART(R,THETA,PHI: Double; var X,Y,Z: Double);
  var RCST : Double;
  begin
    RCST := R*CS(THETA);
    X    := RCST*CS(PHI); Y := RCST*SN(PHI); Z := R*SN(THETA)
  end;

 (*-----------------------------------------------------------------------*)
procedure CROSS(A,B:VECTOR;var C:VECTOR);
  begin
    C[X] := A[Y]*B[Z]-A[Z]*B[Y];
    C[Y] := A[Z]*B[X]-A[X]*B[Z];
    C[Z] := A[X]*B[Y]-A[Y]*B[X];
  end;

(*-----------------------------------------------------------------------*)
function CS(X: Double): Double;
const
  RAD = 0.0174532925199433;
begin
  CS := Cos(X * RAD)
end;

(*-----------------------------------------------------------------------*)
function CUBR(X: Double): Double;
begin
  if (X = 0.0) then
    CUBR := 0.0
  else
    CUBR := EXP(LN(X) / 3.0)
end;

(*-----------------------------------------------------------------------*)
procedure DDD(D, M: integer; S: Double; var DD: Double);
var
  SIGN: Double;
begin
  if ((D < 0) or (M < 0) or (S < 0)) then
    SIGN := -1.0
  else
    SIGN := 1.0;
  DD := SIGN * (abs(D) + abs(M) / 60.0 + abs(S) / 3600.0);
end;

(*-----------------------------------------------------------------------*)

procedure DMS(DD: Double; var D, M: integer; var S: Double);
var
  D1: Double;
begin
  D1 := abs(DD);
  D := Trunc(D1);
  D1 := (D1 - D) * 60.0;
  M := Trunc(D1);
  S := (D1 - M) * 60.0;
  if (DD < 0) then
    if (D <> 0) then
      D := -D
    else if (M <> 0) then
      M := -M
    else
      S := -S;
end;

(*-----------------------------------------------------------------------*)

function DOT(A, B: VECTOR): Double;
begin
  DOT := A[X] * B[X] + A[Y] * B[Y] + A[Z] * B[Z];
end;

(*-----------------------------------------------------------------------*)

procedure LSQFIT(A: LSQMAT; N, M: integer; var S: LSQVEC);
const
  EPS = 1.0E-10; (* machine accuracy *)
var
  I, J, K: integer;
  P, Q, H: Double;
begin
  for J := 1 to M do (* loop over columns 1...M of A *)
    (* eliminate matrix elements A[i,j] with i>j from column j *)
    for I := J + 1 to N do
      if A[I, J] <> 0.0 then
      begin
        (* calculate p, q and new A[j,j]; set A[i,j]=0 *)
        if (abs(A[J, J]) < EPS * abs(A[I, J])) then
        begin
          P := 0.0;
          Q := 1.0;
          A[J, J] := -A[I, J];
          A[I, J] := 0.0;
        end
        else
        begin
          H := Sqrt(A[J, J] * A[J, J] + A[I, J] * A[I, J]);
          if A[J, J] < 0.0 then
            H := -H;
          P := A[J, J] / H;
          Q := -A[I, J] / H;
          A[J, J] := H;
          A[I, J] := 0.0;
        end;
        (* calculate rest of the line *)
        for K := J + 1 to M + 1 do
        begin
          H := P * A[J, K] - Q * A[I, K];
          A[I, K] := Q * A[J, K] + P * A[I, K];
          A[J, K] := H;
        end;
      end;

  (* backsubstitution *)

  for I := M downto 1 do
  begin
    H := A[I, M + 1];
    for K := I + 1 to M do
      H := H - A[I, K] * S[K];
    S[I] := H / A[I, I];
  end;

end; (* LSQFIT *)

(*-----------------------------------------------------------------------*)
function NORM(A: VECTOR): Double;
begin
  NORM := Sqrt(DOT(A, A));
end;

(* ----------------------------------------------------------------------- *)
procedure POLAR(X, Y, Z: Double; var R, THETA, PHI: Double);
var
  RHO: Double;
begin
  RHO := X * X + Y * Y;
  R := Sqrt(RHO + Z * Z);
  PHI := ATN2(Y, X);
  if PHI < 0 then
    PHI := PHI + 360.0;
  RHO := Sqrt(RHO);
  THETA := ATN2(Z, RHO);
end;

(*-----------------------------------------------------------------------*)

procedure QUAD(Y_MINUS, Y_0, Y_PLUS: Double; var XE, YE, ZERO1, ZERO2: Double; var NZ: integer);
var
  A, B, C, DIS, DX: Double;
begin
  NZ := 0;
  A := 0.5 * (Y_MINUS + Y_PLUS) - Y_0;
  B := 0.5 * (Y_PLUS - Y_MINUS);
  C := Y_0;
  XE := -B / (2.0 * A);
  YE := (A * XE + B) * XE + C;
  DIS := B * B - 4.0 * A * C; (* discriminant of y = axx+bx+c *)
  if (DIS >= 0) then (* parabola intersects x-axis *)
  begin
    DX := 0.5 * Sqrt(DIS) / abs(A);
    ZERO1 := XE - DX;
    ZERO2 := XE + DX;
    if (abs(ZERO1) <= +1.0) then
      NZ := NZ + 1;
    if (abs(ZERO2) <= +1.0) then
      NZ := NZ + 1;
    if (ZERO1 < -1.0) then
      ZERO1 := ZERO2;
  end;
end;

(*-----------------------------------------------------------------------*)

function T_EVAL(F: TPolynom; X: Double): Double;
var
  F1, F2, OLD_F1, XX, XX2: Double;
  I: integer;
begin
  if ((X < F.A) or (F.B < X)) then
  begin
    writeln(' T_EVAL : x not within [a,b]');
  end;
  F1 := 0.0;
  F2 := 0.0;
  XX := (2.0 * X - F.A - F.B) / (F.B - F.A);
  XX2 := 2.0 * XX;
  for I := F.M downto 1 do
  begin
    OLD_F1 := F1;
    F1 := XX2 * F1 - F2 + F.C[I];
    F2 := OLD_F1;
  end;
  T_EVAL := XX * F1 - F2 + 0.5 * F.C[0]
end;

(*-----------------------------------------------------------------------*)

procedure T_FIT_LBR(TA, TB: Double; N: integer; var L_POLY, B_POLY, R_POLY: TPolynom);
const
  PI = 3.1415926535898;
  NDIM = 27;
var
  I, J, K: integer;
  FAC, BPA, BMA, PHI: Double;
  T, H, L, B, R: array [0 .. NDIM] of Double;

(*sub*)procedure Positions(T: Double; var LL, BB, RR: Double);
begin
  //
end;

begin
  if (NDIM < 2 * MAX_TP_DEG + 1) then
    writeln(' NDIM too small in T_FIT_LBR');
  if (N > MAX_TP_DEG) then
    writeln(' N too large in T_FIT_LBR');
  L_POLY.M := N;
  B_POLY.M := N;
  R_POLY.M := N;
  L_POLY.A := TA;
  B_POLY.A := TA;
  R_POLY.A := TA;
  L_POLY.B := TB;
  B_POLY.B := TB;
  R_POLY.B := TB;
  BMA := (TB - TA) / 2.0;
  BPA := (TB + TA) / 2.0;
  FAC := 2.0 / (N + 1);
  PHI := PI / (2 * N + 2); (* h(k)=cos(pi*k/N/2) *)
  H[0] := 1.0;
  H[1] := Cos(PHI);
  for I := 2 to (2 * N + 1) do
    H[I] := 2 * H[1] * H[I - 1] - H[I - 2];
  for K := 1 to N + 1 do
    T[K] := H[2 * K - 1] * BMA + BPA; (* subdivison points *)
  for K := 1 to N + 1 do
    Positions(T[K], L[K], B[K], R[K]);
  for K := 2 to N + 1 do (* make L continuous *)
    if (L[K - 1] < L[K]) then
      L[K] := L[K] - 360.0; (* in [-360,+360] ! *)
  for J := 0 to N do (* calculate Chebyshev *)
  begin (* coefficients C(j) *)
    PHI := PI * J / (2 * N + 2);
    H[1] := Cos(PHI);
    for I := 2 to (2 * N + 1) do
      H[I] := 2 * H[1] * H[I - 1] - H[I - 2];
    L_POLY.C[J] := 0.0;
    B_POLY.C[J] := 0.0;
    R_POLY.C[J] := 0.0;
    for K := 1 to N + 1 do
    begin
      L_POLY.C[J] := L_POLY.C[J] + H[2 * K - 1] * L[K];
      B_POLY.C[J] := B_POLY.C[J] + H[2 * K - 1] * B[K];
      R_POLY.C[J] := R_POLY.C[J] + H[2 * K - 1] * R[K];
    end;
    L_POLY.C[J] := L_POLY.C[J] * FAC;
    B_POLY.C[J] := B_POLY.C[J] * FAC;
    R_POLY.C[J] := R_POLY.C[J] * FAC;
  end;
end;

(*-----------------------------------------------------------------------*)

function SN(X: Double): Double;
const
  RAD = 0.0174532925199433;
begin
  SN := Sin(X * RAD)
end;

(*-----------------------------------------------------------------------*)
function TN(X: Double): Double;
const
  RAD = 0.0174532925199433;
var
  XX: Double;
begin
  XX := X * RAD;
  TN := Sin(XX) / Cos(XX);
end;

end.
