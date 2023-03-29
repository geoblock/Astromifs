unit Apc.PrecNut;

interface

uses
  Apc.Mathem;

(*-----------------------------------------------------------------------*)
(* NutEqu: transformation of mean to true coordinates                    *)
(*         (including terms >0.1" according to IAU 1980)                 *)
(*         T = (JD-2451545.0)/36525.0                                    *)
(*-----------------------------------------------------------------------*)
procedure NutEqu(T: double;var X,Y,Z: double);
(*-----------------------------------------------------------------------*)
(* PrecMatEcl: calculates the precession matrix A[i,j] for                  *)
(*          transforming ecliptic coordinates from equinox T1 to T2      *)
(*          ( T=(JD-2451545.0)/36525 )                                   *)
(*-----------------------------------------------------------------------*)
procedure PrecMatEcl(T1,T2: double;var A: Double33);
(*-----------------------------------------------------------------------*)
(* PrecMatEqu: calculates the precession matrix A[i,j] for                  *)
(*          transforming equatorial coordinates from equinox T1 to T2    *)
(*          (T=(JD-2451545.0)/36525 )                                    *)
(*-----------------------------------------------------------------------*)
procedure PrecMatEqu(T1, T2: double; var A: Double33);
(*-----------------------------------------------------------------------*)
(* PN_Matrix: combined precession and nutation matrix for transformation *)
(*            from mean equinox T0 to true equinox T                     *)
(*            T0,T in Julian cent. since J2000; T=(JD-2451545.0)/36525   *)
(*-----------------------------------------------------------------------*)
procedure PN_Matrix(T0, T: double; var A: Double33);
(*-----------------------------------------------------------------------*)
(* PrecArt: calculate change of coordinates due to precession            *)
(*          for given transformation matrix A[i,j]                       *)
(*          (to be used with PrecMatEcl und PrecMatEqu)                  *)
(*-----------------------------------------------------------------------*)
procedure PrecArt(A:Double33; var X,Y,Z: double);

//========================================
implementation
//========================================

(* ----------------------------------------------------------------------- *)
procedure NutEqu(T: double; var X, Y, Z: double);
const
  ARC = 206264.8062; // arcseconds per radian = 3600*180/pi
  P2 = 6.283185307; // 2*pi
var
  LS, D, F, N, EPS: double;
  DPSI, DEPS, C, S: double;
  DX, DY, DZ: double;
  (*
  function Frac(X: double): double;
  // with several compilers it may be necessary to replace TRUNC
  // by LongTrunc or INT if T<-24!
  begin
    Frac := X - Trunc(X)
  end;
  *)
begin
  LS := P2 * Frac(0.993133 + 99.997306 * T); // mean anomaly Sun
  D := P2 * Frac(0.827362 + 1236.853087 * T); // diff. longitude Moon-Sun
  F := P2 * Frac(0.259089 + 1342.227826 * T); // mean argument of latitude
  N := P2 * Frac(0.347346 - 5.372447 * T); // longit. ascending node
  EPS := 0.4090928 - 2.2696E-4 * T; // obliquity of the ecliptic
  DPSI := (-17.200 * Sin(N) - 1.319 * Sin(2 * (F - D + N)) - 0.227 * Sin(2 * (F + N)) + 0.206 *
    Sin(2 * N) + 0.143 * Sin(LS)) / ARC;
  DEPS := (+9.203 * Cos(N) + 0.574 * Cos(2 * (F - D + N)) + 0.098 * Cos(2 * (F + N)) - 0.090 *
    Cos(2 * N)) / ARC;
  C := DPSI * Cos(EPS);
  S := DPSI * Sin(EPS);
  DX := -(C * Y + S * Z);
  DY := (C * X - DEPS * Z);
  DZ := (S * X + DEPS * Y);
  X := X + DX;
  Y := Y + DY;
  Z := Z + DZ;
end;

(* ----------------------------------------------------------------------- *)

procedure PrecMatEcl(T1, T2: double; var A: Double33);
const
  SEC = 3600.0;
var
  DT, PPI, Pi, PA: double;
  C1, S1, C2, S2, C3, S3: double;
begin
  DT := T2 - T1;
  PPI := 174.876383889 + (((3289.4789 + 0.60622 * T1) * T1) + ((-869.8089 - 0.50491 * T1) + 0.03536
    * DT) * DT) / SEC;
  Pi := ((47.0029 - (0.06603 - 0.000598 * T1) * T1) + ((-0.03302 + 0.000598 * T1) + 0.000060 * DT) *
    DT) * DT / SEC;
  PA := ((5029.0966 + (2.22226 - 0.000042 * T1) * T1) + ((1.11113 - 0.000042 * T1) - 0.000006 * DT)
    * DT) * DT / SEC;
  C1 := CS(PPI + PA);
  C2 := CS(Pi);
  C3 := CS(PPI);
  S1 := SN(PPI + PA);
  S2 := SN(Pi);
  S3 := SN(PPI);
  A[1, 1] := +C1 * C3 + S1 * C2 * S3;
  A[1, 2] := +C1 * S3 - S1 * C2 * C3;
  A[1, 3] := -S1 * S2;
  A[2, 1] := +S1 * C3 - C1 * C2 * S3;
  A[2, 2] := +S1 * S3 + C1 * C2 * C3;
  A[2, 3] := +C1 * S2;
  A[3, 1] := +S2 * S3;
  A[3, 2] := -S2 * C3;
  A[3, 3] := +C2;
end;

(* ----------------------------------------------------------------------- *)
procedure PrecMatEqu(T1, T2: double; var A: Double33);
const
  SEC = 3600.0;
var
  DT, ZETA, Z, THETA: double;
  C1, S1, C2, S2, C3, S3: double;
begin
  DT := T2 - T1;
  ZETA := ((2306.2181 + (1.39656 - 0.000139 * T1) * T1) + ((0.30188 - 0.000345 * T1) + 0.017998 *
    DT) * DT) * DT / SEC;
  Z := ZETA + ((0.79280 + 0.000411 * T1) + 0.000205 * DT) * DT * DT / SEC;
  THETA := ((2004.3109 - (0.85330 + 0.000217 * T1) * T1) - ((0.42665 + 0.000217 * T1) + 0.041833 *
    DT) * DT) * DT / SEC;
  C1 := CS(Z);
  C2 := CS(THETA);
  C3 := CS(ZETA);
  S1 := SN(Z);
  S2 := SN(THETA);
  S3 := SN(ZETA);
  A[1, 1] := -S1 * S3 + C1 * C2 * C3;
  A[1, 2] := -S1 * C3 - C1 * C2 * S3;
  A[1, 3] := -C1 * S2;
  A[2, 1] := +C1 * S3 + S1 * C2 * C3;
  A[2, 2] := +C1 * C3 - S1 * C2 * S3;
  A[2, 3] := -S1 * S2;
  A[3, 1] := +S2 * C3;
  A[3, 2] := -S2 * S3;
  A[3, 3] := +C2;
end;

(* ----------------------------------------------------------------------- *)

procedure PN_Matrix(T0, T: double; var A: Double33);
begin
  PrecMatEqu(T0, T, A); // precession matrix T0->T;
  NutEqu(T, A[1, 1], A[2, 1], A[3, 1]); // transform column vectors of
  NutEqu(T, A[1, 2], A[2, 2], A[3, 2]); // matrix A from mean equinox T
  NutEqu(T, A[1, 3], A[2, 3], A[3, 3]); // to true equinox T
end;

(* ----------------------------------------------------------------------- *)

procedure PrecArt(A: Double33; var X, Y, Z: double);
var
  U, V, W: double;
begin
  U := A[1, 1] * X + A[1, 2] * Y + A[1, 3] * Z;
  V := A[2, 1] * X + A[2, 2] * Y + A[2, 3] * Z;
  W := A[3, 1] * X + A[3, 2] * Y + A[3, 3] * Z;
  X := U;
  Y := V;
  Z := W;
end;

end.
