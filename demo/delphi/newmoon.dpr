(* ----------------------------------------------------------------------- *)
(* NewMoon *)
(* Date of New Moon and ecliptic latitude of the Moon *)
(* ----------------------------------------------------------------------- *)

program NewMoon(Input, Output);

{$APPTYPE CONSOLE}

uses
  Apc.Time;

const
  D1 = +1236.853086; (* rate of change dD/dT of the mean elongation *)
  (* of the Moon from the Sun (in revol./century) *)
  D0 = +0.827361; (* mean elongation D of the Moon from the sun *)
  (* for the epoch J2000 in units of 1rev=360deg *)

var
  DAY, MONTH, YEAR, YEAR_CALC: integer;
  HOUR: Double;
  LUNATION_0, LUNATION_I: integer;
  T_NEW_MOON, MJD_NEW_MOON: Double;
  B_MOON: Double;

  (* ----------------------------------------------------------------------- *)
  (* Improve: improves an approximation T for the time of New Moon and *)
  (* finds the ecliptic longitude B of the Moon for that date *)
  (* ( T in julian cent. since J2000, T=(JD-2451545)/36525 ) *)
  (* ----------------------------------------------------------------------- *)
procedure Improve(var T, B: Double);
const
  P2 = 6.283185307; (* 2*pi *)
  ARC = 206264.8062; (* arcsec per radian *)
var
  L, LS, D, F, DLM, DLS, DLAMBDA: Double;
  (* with some compilers it may be necessary to replace TRUNC *)
  (* by LongTrunc or INT if T<-24! *)
  function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1;
    Frac := X
  end;

begin
  (* mean elements L,LS,D,F of the lunar orbit *)
  L := P2 * Frac(0.374897 + 1325.552410 * T); (* mean anomaly of the Moon *)
  LS := P2 * Frac(0.993133 + 99.997361 * T); (* mean anomaly of the Sun *)
  D := P2 * (Frac(0.5 + D0 + D1 * T) - 0.5); (* mean elongation Moon-Sun *)
  F := P2 * Frac(0.259086 + 1342.227825 * T); (* long.Moon-long.asc.node *)
  (* periodic perturbations of the lunar and solar longitude (in") *)
  DLM := +22640 * Sin(L) - 4586 * Sin(L - 2 * D) + 2370 * Sin(2 * D) + 769 * Sin(2 * L) - 668 *
    Sin(LS) - 412 * Sin(2 * F) - 212 * Sin(2 * L - 2 * D) - 206 * Sin(L + LS - 2 * D) + 192 *
    Sin(L + 2 * D) - 165 * Sin(LS - 2 * D) - 125 * Sin(D) - 110 * Sin(L + LS) + 148 * Sin(L - LS) -
    55 * Sin(2 * F - 2 * D);
  DLS := +6893 * Sin(LS) + 72 * Sin(2 * LS);
  (* difference of the true longitudes of Moon and Sun in revolutions *)
  DLAMBDA := D / P2 + (DLM - DLS) / 1296000.0;
  (* correction for the time of newmmon *)
  T := T - DLAMBDA / D1;
  (* ecliptic latitude B of the moon (in deg) *)
  B := (+18520.0 * Sin(F + DLM / ARC) - 526 * Sin(F - 2 * D)) / 3600.0;
end;
(* ----------------------------------------------------------------------- *)

begin (* main program *)

  writeln;
  writeln(' NewMoon: Date of New Moon and ecliptic latitude of the Moon');
  writeln('                       version 93/07/01                     ');
  writeln('         (c) 1993 Thomas Pfleger, Oliver Montenbruck        ');
  writeln;
  write(' Dates of New Moon for the year ');
  readln(YEAR_CALC);
  writeln;
  writeln(' Date ':16, 'UT':7, 'Latitude':12);
  writeln;
  writeln(' h':23, 'o':9);

  LUNATION_0 := Trunc(D1 * (YEAR_CALC - 2000) / 100);
  for LUNATION_I := LUNATION_0 to LUNATION_0 + 13 do
  begin
    T_NEW_MOON := (LUNATION_I - D0) / D1;
    Improve(T_NEW_MOON, B_MOON);
    Improve(T_NEW_MOON, B_MOON);
    MJD_NEW_MOON := 36525.0 * T_NEW_MOON + 51544.5;
    CalDat(MJD_NEW_MOON, DAY, MONTH, YEAR, HOUR);
    if YEAR = YEAR_CALC then
      writeln(YEAR:12, '/', MONTH:2, '/', DAY:2, HOUR:6:1, B_MOON:9:1)
  end;
  writeln;

end.
