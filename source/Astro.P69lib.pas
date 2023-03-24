unit Astro.P69lib;

interface

(*-----------------------------------------------------------------------*)
(* Saturn200: Saturn; ecliptic coordinates L,B,R (in deg and AU)         *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure Saturn200(T: Double; var L, B, R: Double);
(*-----------------------------------------------------------------------*)
(* Uranus200: Uranus; ecliptic coordinates L,B,R (in deg and AU)         *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure Uranus200(T: Double; var L, B, R: Double);
(*-----------------------------------------------------------------------*)
(* Neptune200: Neptune; ecliptic coordinates L,B,R (in deg and AU)       *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure Neptune200(T: Double; var L, B, R: Double);
(*-----------------------------------------------------------------------*)
(* Pluto200: Pluto; ecliptic coordinates L,B,R (in deg and AU)           *)
(*         equinox of date; only valid between 1890 and 2100!!           *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure Pluto200(T: Double; var L, B, R: Double);


implementation

(*-----------------------------------------------------------------------*)
procedure Saturn200(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C6, S6: array [0 .. 11] of Double;
  C, S: array [-6 .. 1] of Double;
  M5, M6, M7, M8: Double;
  U, V, DL, DR, DB: Double;
  I: Integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure ADDTHE(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I6, I, IT: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      ADDTHE(C6[I6], S6[I6], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PERTJUP; (* Kepler terms and perturbations by Jupiter *)
  var
    I: Integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[1] := Cos(M5);
    S[1] := Sin(M5);
    for I := 0 downto -5 do
      ADDTHE(C[I], S[I], C[1], -S[1], C[I - 1], S[I - 1]);
    Term(0, -1, 0, 12.0, -1.4, -13.9, 6.4, 1.2, -1.8);
    Term(0, -2, 0, 0.0, -0.2, -0.9, 1.0, 0.0, -0.1);
    Term(1, 1, 0, 0.9, 0.4, -1.8, 1.9, 0.2, 0.2);
    Term(1, 0, 0, -348.3, 22907.7, -52915.5, -752.2, -3266.5, 8314.4);
    Term(1, 0, 1, -225.2, -146.2, 337.7, -521.3, 79.6, 17.4);
    Term(1, 0, 2, 1.3, -1.4, 3.2, 2.9, 0.1, -0.4);
    Term(1, -1, 0, -1.0, -30.7, 108.6, -815.0, -3.6, -9.3);
    Term(1, -2, 0, -2.0, -2.7, -2.1, -11.9, -0.1, -0.4);
    Term(2, 1, 0, 0.1, 0.2, -1.0, 0.3, 0.0, 0.0);
    Term(2, 0, 0, 44.2, 724.0, -1464.3, -34.7, -188.7, 459.1);
    Term(2, 0, 1, -17.0, -11.3, 18.9, -28.6, 1.0, -3.7);
    Term(2, -1, 0, -3.5, -426.6, -546.5, -26.5, -1.6, -2.7);
    Term(2, -1, 1, 3.5, -2.2, -2.6, -4.3, 0.0, 0.0);
    Term(2, -2, 0, 10.5, -30.9, -130.5, -52.3, -1.9, 0.2);
    Term(2, -3, 0, -0.2, -0.4, -1.2, -0.1, -0.1, 0.0);
    Term(3, 0, 0, 6.5, 30.5, -61.1, 0.4, -11.6, 28.1);
    Term(3, 0, 1, -1.2, -0.7, 1.1, -1.8, -0.2, -0.6);
    Term(3, -1, 0, 29.0, -40.2, 98.2, 45.3, 3.2, -9.4);
    Term(3, -1, 1, 0.6, 0.6, -1.0, 1.3, 0.0, 0.0);
    Term(3, -2, 0, -27.0, -21.1, -68.5, 8.1, -19.8, 5.4);
    Term(3, -2, 1, 0.9, -0.5, -0.4, -2.0, -0.1, -0.8);
    Term(3, -3, 0, -5.4, -4.1, -19.1, 26.2, -0.1, -0.1);
    Term(4, 0, 0, 0.6, 1.4, -3.0, -0.2, -0.6, 1.6);
    Term(4, -1, 0, 1.5, -2.5, 12.4, 4.7, 1.0, -1.1);
    Term(4, -2, 0, -821.9, -9.6, -26.0, 1873.6, -70.5, -4.4);
    Term(4, -2, 1, 4.1, -21.9, -50.3, -9.9, 0.7, -3.0);
    Term(4, -3, 0, -2.0, -4.7, -19.3, 8.2, -0.1, -0.3);
    Term(4, -4, 0, -1.5, 1.3, 6.5, 7.3, 0.0, 0.0);
    Term(5, -2, 0, -2627.6, -1277.3, 117.4, -344.1, -13.8, -4.3);
    Term(5, -2, 1, 63.0, -98.6, 12.7, 6.7, 0.1, -0.2);
    Term(5, -2, 2, 1.7, 1.2, -0.2, 0.3, 0.0, 0.0);
    Term(5, -3, 0, 0.4, -3.6, -11.3, -1.6, 0.0, -0.3);
    Term(5, -4, 0, -1.4, 0.3, 1.5, 6.3, -0.1, 0.0);
    Term(5, -5, 0, 0.3, 0.6, 3.0, -1.7, 0.0, 0.0);
    Term(6, -2, 0, -146.7, -73.7, 166.4, -334.3, -43.6, -46.7);
    Term(6, -2, 1, 5.2, -6.8, 15.1, 11.4, 1.7, -1.0);
    Term(6, -3, 0, 1.5, -2.9, -2.2, -1.3, 0.1, -0.1);
    Term(6, -4, 0, -0.7, -0.2, -0.7, 2.8, 0.0, 0.0);
    Term(6, -5, 0, 0.0, 0.5, 2.5, -0.1, 0.0, 0.0);
    Term(6, -6, 0, 0.3, -0.1, -0.3, -1.2, 0.0, 0.0);
    Term(7, -2, 0, -9.6, -3.9, 9.6, -18.6, -4.7, -5.3);
    Term(7, -2, 1, 0.4, -0.5, 1.0, 0.9, 0.3, -0.1);
    Term(7, -3, 0, 3.0, 5.3, 7.5, -3.5, 0.0, 0.0);
    Term(7, -4, 0, 0.2, 0.4, 1.6, -1.3, 0.0, 0.0);
    Term(7, -5, 0, -0.1, 0.2, 1.0, 0.5, 0.0, 0.0);
    Term(7, -6, 0, 0.2, 0.0, 0.2, -1.0, 0.0, 0.0);
    Term(8, -2, 0, -0.7, -0.2, 0.6, -1.2, -0.4, -0.4);
    Term(8, -3, 0, 0.5, 1.0, -2.0, 1.5, 0.1, 0.2);
    Term(8, -4, 0, 0.4, 1.3, 3.6, -0.9, 0.0, -0.1);
    Term(9, -4, 0, 4.0, -8.7, -19.9, -9.9, 0.2, -0.4);
    Term(9, -4, 1, 0.5, 0.3, 0.8, -1.8, 0.0, 0.0);
    Term(10, -4, 0, 21.3, -16.8, 3.3, 3.3, 0.2, -0.2);
    Term(10, -4, 1, 1.0, 1.7, -0.4, 0.4, 0.0, 0.0);
    Term(11, -4, 0, 1.6, -1.3, 3.0, 3.7, 0.8, -0.2);
  end;

  (*sub*)procedure PERTURA; (* perturbations by Uranus *)
  var
    I: Integer;
  begin
    C[-1] := Cos(M7);
    S[-1] := -Sin(M7);
    for I := -1 downto -4 do
      ADDTHE(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(0, -1, 0, 1.0, 0.7, 0.4, -1.5, 0.1, 0.0);
    Term(0, -2, 0, 0.0, -0.4, -1.1, 0.1, -0.1, -0.1);
    Term(0, -3, 0, -0.9, -1.2, -2.7, 2.1, -0.5, -0.3);
    Term(1, -1, 0, 7.8, -1.5, 2.3, 12.7, 0.0, 0.0);
    Term(1, -2, 0, -1.1, -8.1, 5.2, -0.3, -0.3, -0.3);
    Term(1, -3, 0, -16.4, -21.0, -2.1, 0.0, 0.4, 0.0);
    Term(2, -1, 0, 0.6, -0.1, 0.1, 1.2, 0.1, 0.0);
    Term(2, -2, 0, -4.9, -11.7, 31.5, -13.3, 0.0, -0.2);
    Term(2, -3, 0, 19.1, 10.0, -22.1, 42.1, 0.1, -1.1);
    Term(2, -4, 0, 0.9, -0.1, 0.1, 1.4, 0.0, 0.0);
    Term(3, -2, 0, -0.4, -0.9, 1.7, -0.8, 0.0, -0.3);
    Term(3, -3, 0, 2.3, 0.0, 1.0, 5.7, 0.3, 0.3);
    Term(3, -4, 0, 0.3, -0.7, 2.0, 0.7, 0.0, 0.0);
    Term(3, -5, 0, -0.1, -0.4, 1.1, -0.3, 0.0, 0.0);
  end;

  (*sub*)procedure PERTNEP; (* perturbations by Neptune *)
  begin
    C[-1] := Cos(M8);
    S[-1] := -Sin(M8);
    ADDTHE(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(1, -1, 0, -1.3, -1.2, 2.3, -2.5, 0.0, 0.0);
    Term(1, -2, 0, 1.0, -0.1, 0.1, 1.4, 0.0, 0.0);
    Term(2, -2, 0, 1.1, -0.1, 0.2, 3.3, 0.0, 0.0);
  end;

  (*sub*)procedure PERTJUR; (* perturbations by Jupiter and Uranus *)
  var
    PHI, X, Y: Double;
  begin
    PHI := (-2 * M5 + 5 * M6 - 3 * M7);
    X := Cos(PHI);
    Y := Sin(PHI);
    DL := DL - 0.8 * X - 0.1 * Y;
    DR := DR - 0.2 * X + 1.8 * Y;
    DB := DB + 0.3 * X + 0.5 * Y;
    ADDTHE(X, Y, C6[1], S6[1], X, Y);
    DL := DL + (+2.4 - 0.7 * T) * X + (27.8 - 0.4 * T) * Y;
    DR := DR + 2.1 * X - 0.2 * Y;
    ADDTHE(X, Y, C6[1], S6[1], X, Y);
    DL := DL + 0.1 * X + 1.6 * Y;
    DR := DR - 3.6 * X + 0.3 * Y;
    DB := DB - 0.2 * X + 0.6 * Y;
  end;

begin (* Saturn200 *)

  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M5 := P2 * Frac(0.0565314 + 8.4302963 * T);
  M6 := P2 * Frac(0.8829867 + 3.3947688 * T);
  M7 := P2 * Frac(0.3969537 + 1.1902586 * T);
  M8 := P2 * Frac(0.7208473 + 0.6068623 * T);
  C6[0] := 1.0;
  S6[0] := 0.0;
  C6[1] := Cos(M6);
  S6[1] := Sin(M6);
  for I := 2 to 11 do
    ADDTHE(C6[I - 1], S6[I - 1], C6[1], S6[1], C6[I], S6[I]);
  PERTJUP;
  PERTURA;
  PERTNEP;
  PERTJUR;
  L := 360.0 * Frac(0.2561136 + M6 / P2 + ((5018.6 + T * 1.9) * T + DL) / 1296.0E3);
  R := 9.557584 - 0.000186 * T + DR * 1.0E-5;
  B := (175.1 - 10.2 * T + DB) / 3600.0;

end; (* Saturn200 *)

(* ----------------------------------------------------------------------- *)

procedure Uranus200(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C7, S7: array [-2 .. 7] of Double;
  C, S: array [-8 .. 0] of Double;
  M5, M6, M7, M8: Double;
  U, V, DL, DR, DB: Double;
  I: Integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure ADDTHE(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I7, I, IT: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      ADDTHE(C7[I7], S7[I7], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PERTJUP; (* perturbations by Jupiter *)
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    ADDTHE(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(-1, -1, 0, 0.0, 0.0, -0.1, 1.7, -0.1, 0.0);
    Term(0, -1, 0, 0.5, -1.2, 18.9, 9.1, -0.9, 0.1);
    Term(1, -1, 0, -21.2, 48.7, -455.5, -198.8, 0.0, 0.0);
    Term(1, -2, 0, -0.5, 1.2, -10.9, -4.8, 0.0, 0.0);
    Term(2, -1, 0, -1.3, 3.2, -23.2, -11.1, 0.3, 0.1);
    Term(2, -2, 0, -0.2, 0.2, 1.1, 1.5, 0.0, 0.0);
    Term(3, -1, 0, 0.0, 0.2, -1.8, 0.4, 0.0, 0.0);
  end;

  (*sub*)procedure PERTSAT; (* perturbations by Saturn *)
  var
    I: Integer;
  begin
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    for I := -1 downto -3 do
      ADDTHE(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(0, -1, 0, 1.4, -0.5, -6.4, 9.0, -0.4, -0.8);
    Term(1, -1, 0, -18.6, -12.6, 36.7, -336.8, 1.0, 0.3);
    Term(1, -2, 0, -0.7, -0.3, 0.5, -7.5, 0.1, 0.0);
    Term(2, -1, 0, 20.0, -141.6, -587.1, -107.0, 3.1, -0.8);
    Term(2, -1, 1, 1.0, 1.4, 5.8, -4.0, 0.0, 0.0);
    Term(2, -2, 0, 1.6, -3.8, -35.6, -16.0, 0.0, 0.0);
    Term(3, -1, 0, 75.3, -100.9, 128.9, 77.5, -0.8, 0.1);
    Term(3, -1, 1, 0.2, 1.8, -1.9, 0.3, 0.0, 0.0);
    Term(3, -2, 0, 2.3, -1.3, -9.5, -17.9, 0.0, 0.1);
    Term(3, -3, 0, -0.7, -0.5, -4.9, 6.8, 0.0, 0.0);
    Term(4, -1, 0, 3.4, -5.0, 21.6, 14.3, -0.8, -0.5);
    Term(4, -2, 0, 1.9, 0.1, 1.2, -12.1, 0.0, 0.0);
    Term(4, -3, 0, -0.1, -0.4, -3.9, 1.2, 0.0, 0.0);
    Term(4, -4, 0, -0.2, 0.1, 1.6, 1.8, 0.0, 0.0);
    Term(5, -1, 0, 0.2, -0.3, 1.0, 0.6, -0.1, 0.0);
    Term(5, -2, 0, -2.2, -2.2, -7.7, 8.5, 0.0, 0.0);
    Term(5, -3, 0, 0.1, -0.2, -1.4, -0.4, 0.0, 0.0);
    Term(5, -4, 0, -0.1, 0.0, 0.1, 1.2, 0.0, 0.0);
    Term(6, -2, 0, -0.2, -0.6, 1.4, -0.7, 0.0, 0.0);
  end;

  (*sub*)procedure PERTNEP; (* Kepler terms and perturbations by Neptune *)
  var
    I: Integer;
  begin
    C[-1] := Cos(M8);
    S[-1] := -Sin(M8);
    for I := -1 downto -7 do
      ADDTHE(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, 0, 0, -78.1, 19518.1, -90718.2, -334.7, 2759.5, -311.9);
    Term(1, 0, 1, -81.6, 107.7, -497.4, -379.5, -2.8, -43.7);
    Term(1, 0, 2, -6.6, -3.1, 14.4, -30.6, -0.4, -0.5);
    Term(1, 0, 3, 0.0, -0.5, 2.4, 0.0, 0.0, 0.0);
    Term(2, 0, 0, -2.4, 586.1, -2145.2, -15.3, 130.6, -14.3);
    Term(2, 0, 1, -4.5, 6.6, -24.2, -17.8, 0.7, -1.6);
    Term(2, 0, 2, -0.4, 0.0, 0.1, -1.4, 0.0, 0.0);
    Term(3, 0, 0, 0.0, 24.5, -76.2, -0.6, 7.0, -0.7);
    Term(3, 0, 1, -0.2, 0.4, -1.4, -0.8, 0.1, -0.1);
    Term(4, 0, 0, 0.0, 1.1, -3.0, 0.1, 0.4, 0.0);
    Term(-1, -1, 0, -0.2, 0.2, 0.7, 0.7, -0.1, 0.0);
    Term(0, -1, 0, -2.8, 2.5, 8.7, 10.5, -0.4, -0.1);
    Term(1, -1, 0, -28.4, 20.3, -51.4, -72.0, 0.0, 0.0);
    Term(1, -2, 0, -0.6, -0.1, 4.2, -14.6, 0.2, 0.4);
    Term(1, -3, 0, 0.2, 0.5, 3.4, -1.6, -0.1, 0.1);
    Term(2, -1, 0, -1.8, 1.3, -5.5, -7.7, 0.0, 0.3);
    Term(2, -2, 0, 29.4, 10.2, -29.0, 83.2, 0.0, 0.0);
    Term(2, -3, 0, 8.8, 17.8, -41.9, 21.5, -0.1, -0.3);
    Term(2, -4, 0, 0.0, 0.1, -2.1, -0.9, 0.1, 0.0);
    Term(3, -2, 0, 1.5, 0.5, -1.7, 5.1, 0.1, -0.2);
    Term(3, -3, 0, 4.4, 14.6, -84.3, 25.2, 0.1, -0.1);
    Term(3, -4, 0, 2.4, -4.5, 12.0, 6.2, 0.0, 0.0);
    Term(3, -5, 0, 2.9, -0.9, 2.1, 6.2, 0.0, 0.0);
    Term(4, -3, 0, 0.3, 1.0, -4.0, 1.1, 0.1, -0.1);
    Term(4, -4, 0, 2.1, -2.7, 17.9, 14.0, 0.0, 0.0);
    Term(4, -5, 0, 3.0, -0.4, 2.3, 17.6, -0.1, -0.1);
    Term(4, -6, 0, -0.6, -0.5, 1.1, -1.6, 0.0, 0.0);
    Term(5, -4, 0, 0.2, -0.2, 1.0, 0.8, 0.0, 0.0);
    Term(5, -5, 0, -0.9, -0.1, 0.6, -7.1, 0.0, 0.0);
    Term(5, -6, 0, -0.5, -0.6, 3.8, -3.6, 0.0, 0.0);
    Term(5, -7, 0, 0.0, -0.5, 3.0, 0.1, 0.0, 0.0);
    Term(6, -6, 0, 0.2, 0.3, -2.7, 1.6, 0.0, 0.0);
    Term(6, -7, 0, -0.1, 0.2, -2.0, -0.4, 0.0, 0.0);
    Term(7, -7, 0, 0.1, -0.2, 1.3, 0.5, 0.0, 0.0);
    Term(7, -8, 0, 0.1, 0.0, 0.4, 0.9, 0.0, 0.0);
  end;

  (*sub*)procedure PERTJSU; (* perturbations by Jupiter and Saturn *)
  var
    I: Integer;
  begin
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    C[-4] := Cos(-4 * M6 + 2 * M5);
    S[-4] := Sin(-4 * M6 + 2 * M5);
    for I := -4 downto -5 do
      ADDTHE(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(-2, -4, 0, -0.7, 0.4, -1.5, -2.5, 0.0, 0.0);
    Term(-1, -4, 0, -0.1, -0.1, -2.2, 1.0, 0.0, 0.0);
    Term(1, -5, 0, 0.1, -0.4, 1.4, 0.2, 0.0, 0.0);
    Term(1, -6, 0, 0.4, 0.5, -0.8, -0.8, 0.0, 0.0);
    Term(2, -6, 0, 5.7, 6.3, 28.5, -25.5, 0.0, 0.0);
    Term(2, -6, 1, 0.1, -0.2, -1.1, -0.6, 0.0, 0.0);
    Term(3, -6, 0, -1.4, 29.2, -11.4, 1.1, 0.0, 0.0);
    Term(3, -6, 1, 0.8, -0.4, 0.2, 0.3, 0.0, 0.0);
    Term(4, -6, 0, 0.0, 1.3, -6.0, -0.1, 0.0, 0.0);
  end;

begin (* Uranus200 *)

  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M5 := P2 * Frac(0.0564472 + 8.4302889 * T);
  M6 := P2 * Frac(0.8829611 + 3.3947583 * T);
  M7 := P2 * Frac(0.3967117 + 1.1902849 * T);
  M8 := P2 * Frac(0.7216833 + 0.6068528 * T);
  C7[0] := 1.0;
  S7[0] := 0.0;
  C7[1] := Cos(M7);
  S7[1] := Sin(M7);
  for I := 2 to 7 do
    ADDTHE(C7[I - 1], S7[I - 1], C7[1], S7[1], C7[I], S7[I]);
  for I := 1 to 2 do
  begin
    C7[-I] := C7[I];
    S7[-I] := -S7[I];
  end;
  PERTJUP;
  PERTSAT;
  PERTNEP;
  PERTJSU;
  L := 360.0 * Frac(0.4734843 + M7 / P2 + ((5082.3 + 34.2 * T) * T + DL) / 1296.0E3);
  R := 19.211991 + (-0.000333 - 0.000005 * T) * T + DR * 1.0E-5;
  B := (-130.61 + (-0.54 + 0.04 * T) * T + DB) / 3600.0;

end; (* Uranus200 *)

(* ----------------------------------------------------------------------- *)

procedure Neptune200(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C8, S8: array [0 .. 6] of Double;
  C, S: array [-6 .. 0] of Double;
  M5, M6, M7, M8: Double;
  U, V, DL, DR, DB: Double;
  I: Integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure ADDTHE(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I1, I, IT: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      ADDTHE(C8[I1], S8[I1], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PERTJUP; (* perturbations by Jupiter *)
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    ADDTHE(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(0, -1, 0, 0.1, 0.1, -3.0, 1.8, -0.3, -0.3);
    Term(1, 0, 0, 0.0, 0.0, -15.9, 9.0, 0.0, 0.0);
    Term(1, -1, 0, -17.6, -29.3, 416.1, -250.0, 0.0, 0.0);
    Term(1, -2, 0, -0.4, -0.7, 10.4, -6.2, 0.0, 0.0);
    Term(2, -1, 0, -0.2, -0.4, 2.4, -1.4, 0.4, -0.3);
  end;

  (*sub*)procedure PERTSAT; (* perturbations by Saturn *)
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    ADDTHE(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(0, -1, 0, -0.1, 0.0, 0.2, -1.8, -0.1, -0.5);
    Term(1, 0, 0, 0.0, 0.0, -8.3, -10.4, 0.0, 0.0);
    Term(1, -1, 0, 13.6, -12.7, 187.5, 201.1, 0.0, 0.0);
    Term(1, -2, 0, 0.4, -0.4, 4.5, 4.5, 0.0, 0.0);
    Term(2, -1, 0, 0.4, -0.1, 1.7, -3.2, 0.2, 0.2);
    Term(2, -2, 0, -0.1, 0.0, -0.2, 2.7, 0.0, 0.0);
  end;

  (*sub*)procedure PERTURA; (* Kepler terms and perturbations by Uranus *)
  var
    I: Integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M7);
    S[-1] := -Sin(M7);
    for I := -1 downto -5 do
      ADDTHE(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, 0, 0, 32.3, 3549.5, -25880.2, 235.8, -6360.5, 374.0);
    Term(1, 0, 1, 31.2, 34.4, -251.4, 227.4, 34.9, 29.3);
    Term(1, 0, 2, -1.4, 3.9, -28.6, -10.1, 0.0, -0.9);
    Term(2, 0, 0, 6.1, 68.0, -111.4, 2.0, -54.7, 3.7);
    Term(2, 0, 1, 0.8, -0.2, -2.1, 2.0, -0.2, 0.8);
    Term(3, 0, 0, 0.1, 1.0, -0.7, 0.0, -0.8, 0.1);
    Term(0, -1, 0, -0.1, -0.3, -3.6, 0.0, 0.0, 0.0);
    Term(1, 0, 0, 0.0, 0.0, 5.5, -6.9, 0.1, 0.0);
    Term(1, -1, 0, -2.2, -1.6, -116.3, 163.6, 0.0, -0.1);
    Term(1, -2, 0, 0.2, 0.1, -1.2, 0.4, 0.0, -0.1);
    Term(2, -1, 0, 4.2, -1.1, -4.4, -34.6, -0.2, 0.1);
    Term(2, -2, 0, 8.6, -2.9, -33.4, -97.0, 0.2, 0.1);
    Term(3, -1, 0, 0.1, -0.2, 2.1, -1.2, 0.0, 0.1);
    Term(3, -2, 0, -4.6, 9.3, 38.2, 19.8, 0.1, 0.1);
    Term(3, -3, 0, -0.5, 1.7, 23.5, 7.0, 0.0, 0.0);
    Term(4, -2, 0, 0.2, 0.8, 3.3, -1.5, -0.2, -0.1);
    Term(4, -3, 0, 0.9, 1.7, 17.9, -9.1, -0.1, 0.0);
    Term(4, -4, 0, -0.4, -0.4, -6.2, 4.8, 0.0, 0.0);
    Term(5, -3, 0, -1.6, -0.5, -2.2, 7.0, 0.0, 0.0);
    Term(5, -4, 0, -0.4, -0.1, -0.7, 5.5, 0.0, 0.0);
    Term(5, -5, 0, 0.2, 0.0, 0.0, -3.5, 0.0, 0.0);
    Term(6, -4, 0, -0.3, 0.2, 2.1, 2.7, 0.0, 0.0);
    Term(6, -5, 0, 0.1, -0.1, -1.4, -1.4, 0.0, 0.0);
    Term(6, -6, 0, -0.1, 0.1, 1.4, 0.7, 0.0, 0.0);
  end;

begin (* Neptune200 *)

  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M5 := P2 * Frac(0.0563867 + 8.4298907 * T);
  M6 := P2 * Frac(0.8825086 + 3.3957748 * T);
  M7 := P2 * Frac(0.3965358 + 1.1902851 * T);
  M8 := P2 * Frac(0.7214906 + 0.6068526 * T);
  C8[0] := 1.0;
  S8[0] := 0.0;
  C8[1] := Cos(M8);
  S8[1] := Sin(M8);
  for I := 2 to 6 do
    ADDTHE(C8[I - 1], S8[I - 1], C8[1], S8[1], C8[I], S8[I]);
  PERTJUP;
  PERTSAT;
  PERTURA;
  L := 360.0 * Frac(0.1254046 + M8 / P2 + ((4982.8 - 21.3 * T) * T + DL) / 1296.0E3);
  R := 30.072984 + (0.001234 + 0.000003 * T) * T + DR * 1.0E-5;
  B := (54.77 + (0.26 + 0.06 * T) * T + DB) / 3600.0;

end; (* Neptune200 *)

(* ----------------------------------------------------------------------- *)
procedure Pluto200(T: Double; var L, B, R: Double);

const
  P2 = 6.283185307;
var
  C9, S9: array [0 .. 6] of Double;
  C, S: array [-3 .. 2] of Double;
  M5, M6, M9: Double;
  DL, DR, DB: Double;
  I: Integer;

  (* sub *) function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

(* sub *) procedure ADDTHE(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

(* sub *) procedure Term(I9, I: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  var
    U, V: Double;
  begin
    ADDTHE(C9[I9], S9[I9], C[I], S[I], U, V);
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

(* sub *) procedure PERTJUP; (* Kepler terms and perturbations by Jupiter *)
  var
    I: Integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[1] := Cos(M5);
    S[1] := Sin(M5);
    for I := 0 downto -1 do
      ADDTHE(C[I], S[I], C[1], -S[1], C[I - 1], S[I - 1]);
    ADDTHE(C[1], S[1], C[1], S[1], C[2], S[2]);
    Term(1, 0, 0.06, 100924.08, -960396.0, 15965.1, 51987.68, -24288.76);
    Term(2, 0, 3274.74, 17835.12, -118252.2, 3632.4, 12687.49, -6049.72);
    Term(3, 0, 1543.52, 4631.99, -21446.6, 1167.0, 3504.00, -1853.10);
    Term(4, 0, 688.99, 1227.08, -4823.4, 213.5, 1048.19, -648.26);
    Term(5, 0, 242.27, 415.93, -1075.4, 140.6, 302.33, -209.76);
    Term(6, 0, 138.41, 110.91, -308.8, -55.3, 109.52, -93.82);
    Term(3, -1, -0.99, 5.06, -25.6, 19.8, 1.26, -1.96);
    Term(2, -1, 7.15, 5.61, -96.7, 57.2, 1.64, -2.16);
    Term(1, -1, 10.79, 23.13, -390.4, 236.4, -0.33, 0.86);
    Term(0, 1, -0.23, 4.43, 102.8, 63.2, 3.15, 0.34);
    Term(1, 1, -1.10, -0.92, 11.8, -2.3, 0.43, 0.14);
    Term(2, 1, 0.62, 0.84, 2.3, 0.7, 0.05, -0.04);
    Term(3, 1, -0.38, -0.45, 1.2, -0.8, 0.04, 0.05);
    Term(4, 1, 0.17, 0.25, 0.0, 0.2, -0.01, -0.01);
    Term(3, -2, 0.06, 0.07, -0.6, 0.3, 0.03, -0.03);
    Term(2, -2, 0.13, 0.20, -2.2, 1.5, 0.03, -0.07);
    Term(1, -2, 0.32, 0.49, -9.4, 5.7, -0.01, 0.03);
    Term(0, -2, -0.04, -0.07, 2.6, -1.5, 0.07, -0.02);
  end;

(* sub *) procedure PERTSAT; (* perturbations by Saturn *)
  var
    I: Integer;
  begin
    C[1] := Cos(M6);
    S[1] := Sin(M6);
    for I := 0 downto -1 do
      ADDTHE(C[I], S[I], C[1], -S[1], C[I - 1], S[I - 1]);
    Term(1, -1, -29.47, 75.97, -106.4, -204.9, -40.71, -17.55);
    Term(0, 1, -13.88, 18.20, 42.6, -46.1, 1.13, 0.43);
    Term(1, 1, 5.81, -23.48, 15.0, -6.8, -7.48, 3.07);
    Term(2, 1, -10.27, 14.16, -7.9, 0.4, 2.43, -0.09);
    Term(3, 1, 6.86, -10.66, 7.3, -0.3, -2.25, 0.69);
    Term(2, -2, 4.32, 2.00, 0.0, -2.2, -0.24, 0.12);
    Term(1, -2, -5.04, -0.83, -9.2, -3.1, 0.79, -0.24);
    Term(0, -2, 4.25, 2.48, -5.9, -3.3, 0.58, 0.02);
  end;

(* sub *) procedure PERTJUS; (* perturbations by Jupiter and Saturn *)
  var
    PHI, X, Y: Double;
  begin
    PHI := (M5 - M6);
    X := Cos(PHI);
    Y := Sin(PHI);
    DL := DL - 9.11 * X + 0.12 * Y;
    DR := DR - 3.4 * X - 3.3 * Y;
    DB := DB + 0.81 * X + 0.78 * Y;
    ADDTHE(X, Y, C9[1], S9[1], X, Y);
    DL := DL + 5.92 * X + 0.25 * Y;
    DR := DR + 2.3 * X - 3.8 * Y;
    DB := DB - 0.67 * X - 0.51 * Y;
  end;

(* sub *) procedure PREC(T: Double; var L, B: Double); (* precess. 1950->equinox of date *)
  const
    DEG = 57.2957795;
  var
    D, PPI, PI, P, C1, S1, C2, S2, C3, S3, X, Y, Z: Double;
  begin
    D := T + 0.5;
    L := L / DEG;
    B := B / DEG;
    PPI := 3.044;
    PI := 2.28E-4 * D;
    P := (0.0243764 + 5.39E-6 * D) * D;
    C1 := Cos(PI);
    C2 := Cos(B);
    C3 := Cos(PPI - L);
    S1 := Sin(PI);
    S2 := Sin(B);
    S3 := Sin(PPI - L);
    X := C2 * C3;
    Y := C1 * C2 * S3 - S1 * S2;
    Z := S1 * C2 * S3 + C1 * S2;
    B := DEG * ArcTan(Z / Sqrt((1.0 - Z) * (1.0 + Z)));
    if (X > 0) then
      L := 360.0 * Frac((PPI + P - ArcTan(Y / X)) / P2)
    else
      L := 360.0 * Frac((PPI + P - ArcTan(Y / X)) / P2 + 0.5);
  end;

begin (* Pluto200 *)

  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M5 := P2 * Frac(0.0565314 + 8.4302963 * T);
  M6 := P2 * Frac(0.8829867 + 3.3947688 * T);
  M9 := P2 * Frac(0.0385795 + 0.4026667 * T);
  C9[0] := 1.0;
  S9[0] := 0.0;
  C9[1] := Cos(M9);
  S9[1] := Sin(M9);
  for I := 2 to 6 do
    ADDTHE(C9[I - 1], S9[I - 1], C9[1], S9[1], C9[I], S9[I]);
  PERTJUP;
  PERTSAT;
  PERTJUS;
  L := 360.0 * Frac(0.6232469 + M9 / P2 + DL / 1296.0E3);
  R := 40.7247248 + DR * 1.0E-5;
  B := -3.909434 + DB / 3600.0;
  PREC(T, L, B);

end; (* Pluto200 *)




end.