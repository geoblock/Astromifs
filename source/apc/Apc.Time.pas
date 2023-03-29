unit Apc.Time;

//
//    Functions and classes handling time
//

interface

uses
  System.Math,
  Apc.Mathem;

(*-----------------------------------------------------------------------*)
(* ETminUT: Difference of Ephemeris Time and Universal Time              *)
(*          (polynomial approximation valid from 1900 to 1995)           *)
(*          T:     time in Julian centuries since J2000                  *)
(*                 ( = (JD-2451545.0)/36525.0 )                          *)
(*          DTSEC: DT=ET-UT in sec (only if VALID=True)                  *)
(*          VALID: True for times between 1900 and 1995, False otherwise *)
(*-----------------------------------------------------------------------*)
procedure ETminUT(T: Double; var DTSEC: Double; var VALID: Boolean);

//------------------------------------------------------------------------
//
// GMST: Greenwich mean sidereal time
//       Time as Modified Julian Date
//
// result:   GMST in [rad]
//
//------------------------------------------------------------------------
function GMST(MJD: double): Double;

//-----------------------------------------------------------------------
// LMST: local mean sidereal time
//-----------------------------------------------------------------------
function LMST(MJD, LAMBDA: Double): Double;

(*-----------------------------------------------------------------------*)
(* MJD: Modified Julian Date                                             *)
(*      The routine is valid for any date since 4713 BC.                 *)
(*      Julian calendar is used up to 1582 October 4,                    *)
(*      Gregorian calendar is used from 1582 October 15 onwards.         *)
(*-----------------------------------------------------------------------*)
function MJD(DAY, MONTH, YEAR: integer; HOUR: Double): Double;

(*----------------------------------------------------------------------*)
(* CalDat: Finds the civil calendar date for a given value              *)
(*         of the Modified Julian Date (MJD).                           *)
(*         Julian calendar is used up to 1582 October 4,                *)
(*         Gregorian calendar is used from 1582 October 15 onwards.     *)
(*----------------------------------------------------------------------*)
procedure CalDat(MJD: Double; var DAY, MONTH, YEAR: integer; var HOUR: Double);

//==============================================================
implementation
//==============================================================

function LongTrunc(X: Double): Double;
begin
  LongTrunc := Int(X);
end;

//----------------------------------------------------------------------
procedure CalDat(MJD: Double; var DAY, MONTH, YEAR: integer; var HOUR: Double);
var
  B, D, F: integer;
  JD, JD0, C, E: Double;
begin
  JD := MJD + 2400000.5;
  JD0 := Int(JD + 0.5);  // Delphi
  (* JD0 := Trunc(JD + 0.5); *)  // Standard Pascal
  (* JD0 := LongTrunc(JD + 0.5); *)
  if (JD0 < 2299161.0) // calendar:
  then
  begin
    C := JD0 + 1524.0
  end // -> Julian
  else
  begin // -> Gregorian
    B := Trunc((JD0 - 1867216.25) / 36524.25);
    C := JD0 + (B - Trunc(B / 4)) + 1525.0
  end;
  D := Trunc((C - 122.1) / 365.25);
  E := 365.0 * D + Trunc(D / 4);
  F := Trunc((C - E) / 30.6001);
  DAY := Trunc(C - E + 0.5) - Trunc(30.6001 * F);
  MONTH := F - 1 - 12 * Trunc(F / 14);
  YEAR := D - 4715 - Trunc((7 + MONTH) / 10);
  HOUR := 24.0 * (JD + 0.5 - JD0);
end;  // CalDat

//-----------------------------------------------------------------------

procedure ETminUT(T: Double; var DTSEC: Double; var VALID: Boolean);
begin
  VALID := ((-1.0 <= T) and (T <= -0.05));
  if (VALID) then
  begin
    DTSEC := ((((-449.50 * T - 783.42) * T - 387.70) * T) + 13.34) * T + 62.14;
  end;
end; // ETminUT


//------------------------------------------------------------------------
function GMST(MJD: Double): Double;

const
  Secs = 86400; // Seconds per day
var
  MJD_0, UT, T_0, T, GMST: Double;
  Modulo: Integer;

begin
  MJD_0 := floor(MJD);
  UT := Secs * (MJD - MJD_0); // [s]
  T_0 := (MJD_0 - 51544.5) / 36525.0;
  T := (MJD - 51544.5) / 36525.0;

  GMST := 24110.54841 + 8640184.812866 * T_0 + 1.0027379093 * UT + (0.093104 - 6.2E-6 * T) * T * T;
  // [sec]
  Modulo := Round(GMST) mod Secs;
  // c:-> return (pi2/Secs)*Modulo(gmst,Secs);   // [Rad]
  result := (2 * pi / Secs) * Modulo; // Modulo(gmst,Secs) in [Rad]
end;

(* ----------------------------------------------------------------------- *)
function LMST(MJD, LAMBDA: Double): Double;
var
  MJD0, T, UT, GMST: Double;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1;
    Frac := X
  end;

begin
  (* MJD0:=Trunc(MJD); *)       (* Standard Pascal *)
  MJD0 := Int(MJD); (* Delphi *)
  (* MJD0:=LongTrunc(MJD); *)  (* ST Pascal plus *)
  UT := (MJD - MJD0) * 24;
  T := (MJD0 - 51544.5) / 36525.0;
  GMST := 6.697374558 + 1.0027379093 * UT + (8640184.812866 + (0.093104 - 6.2E-6 * T) * T) *
    T / 3600.0;
  LMST := 24.0 * Frac((GMST - LAMBDA / 15.0) / 24.0);
end;

(* ----------------------------------------------------------------------- *)
function MJD(DAY, MONTH, YEAR: integer; HOUR: Double): Double;
var
  A: Double;
  B: integer;
begin
  A := 10000.0 * YEAR + 100.0 * MONTH + DAY;
  if (MONTH <= 2) then
  begin
    MONTH := MONTH + 12;
    YEAR := YEAR - 1
  end;
  if (A <= 15821004.1) then
    B := -2 + Trunc((YEAR + 4716) / 4) - 1179
  else
    B := Trunc(YEAR / 400) - Trunc(YEAR / 100) + Trunc(YEAR / 4);
  A := 365.0 * YEAR - 679004.0;
  MJD := A + B + Trunc(30.6001 * (MONTH + 1)) + DAY + HOUR / 24.0;
end;  // MJD

end.
