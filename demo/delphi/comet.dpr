(* ----------------------------------------------------------------------- *)
(* Comet *)
(* calculation of unperturbed ephemeris with arbitrary eccentricity *)
(* for comets and minor planets *)
(* ----------------------------------------------------------------------- *)

program Comet(Input, Output, COMINP);

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric,
  Apc.Sun,
  Apc.Kepler,
  Apc.Time;

var
  DAY, MONTH, YEAR, NLINE: integer;
  D, HOUR, T0, Q, ECC, TEQX0, FAC: Double;
  MODJD, T, T1, DT, T2, TEQX: Double;
  X, Y, Z, VX, VY, VZ, XS, YS, ZS: Double;
  L, B, R, LS, BS, RS, Ra, Dec, DELTA, DELTA0: Double;
  PQR, A, ASI: Double33;
  COMINP: TEXT;

  (* ----------------------------------------------------------------------- *)
  (* GetElm: reads orbital elements from file COMINP *)
  (* ----------------------------------------------------------------------- *)
procedure GetElm(var T0, Q, ECC: Double; var PQR: Double33; var TEQX0: Double);
var
  INC, LAN, AOP: Double;
begin
  writeln;
  writeln(' Comet: ephemeris calculation for comets and minor planets');
  writeln('                     version 93/07/01                     ');
  writeln('        (C) 1993 Thomas Pfleger, Oliver Montenbruck       ');
  writeln;
  writeln(' Orbital elements from file COMINP: ');
  writeln;

  (* open file for reading *)
  (* Reset(COMINP); *)                              (* Standard Pascal *)
  Assign(COMINP, 'COMINP.DAT');
  Reset(COMINP);  

  (* display orbital elements *)
  Readln(COMINP, YEAR, MONTH, D);
  writeln('  perihelion time (y m d) ', YEAR:7, MONTH:3, D:6:2);
  Readln(COMINP, Q);
  writeln('  perihelion distance (q) ', Q:14:7, ' AU  ');
  Readln(COMINP, ECC);
  writeln('  eccentricity (e)        ', ECC:14:7);
  Readln(COMINP, INC);
  writeln('  inclination (i)         ', INC:12:5, ' deg');
  Readln(COMINP, LAN);
  writeln('  long. of ascending node ', LAN:12:5, ' deg');
  Readln(COMINP, AOP);
  writeln('  argument of perihelion  ', AOP:12:5, ' deg');
  Readln(COMINP, TEQX0);
  writeln('  equinox                 ', TEQX0:9:2);

  DAY := Trunc(D);
  HOUR := 24.0 * (D - DAY);
  T0 := (MJD(DAY, MONTH, YEAR, HOUR) - 51544.5) / 36525.0;
  TEQX0 := (TEQX0 - 2000.0) / 100.0;
  GaussVec(LAN, INC, AOP, PQR);
end;

(* ----------------------------------------------------------------------- *)
(* GetEph: reads desired period of time and equinox of the ephemeris *)
(* ----------------------------------------------------------------------- *)
procedure GetEph(var T1, DT, T2, TEQX: Double);

var
  YEAR, MONTH, DAY: integer;
  EQX, HOUR, JD: Double;

begin

  writeln;
  writeln(' Begin and end of the ephemeris: ');
  writeln;
  write('  first date (yyyy mm dd hh.hhh)            ... ');
  Readln(YEAR, MONTH, DAY, HOUR);
  T1 := (MJD(DAY, MONTH, YEAR, HOUR) - 51544.5) / 36525.0;
  write('  final date (yyyy mm dd hh.hhh)            ... ');
  Readln(YEAR, MONTH, DAY, HOUR);
  T2 := (MJD(DAY, MONTH, YEAR, HOUR) - 51544.5) / 36525.0;
  write('  step size (dd hh.hh)                      ... ');
  Readln(DAY, HOUR);
  DT := (DAY + HOUR / 24.0) / 36525.0;
  writeln;
  write(' Desired equinox of the ephemeris (yyyy.y)  ... ');
  Readln(EQX);
  TEQX := (EQX - 2000.0) / 100.0;

  writeln;
  writeln;
  writeln('    Date      ET   Sun     l      b     r', '        Ra          Dec      Distance ');
  writeln(' ':45, '   h  m  s      o  ''  "     (AU) ');

end;

(* ----------------------------------------------------------------------- *)
(* WrtLBR: write L and B in deg,min,sec and R *)
(* ----------------------------------------------------------------------- *)
procedure WrtLBR(L, B, R: Double);
var
  H, M: integer;
  S: Double;
begin
  DMS(L, H, M, S);
  write(H:5, M:3, S:5:1);
  DMS(B, H, M, S);
  writeln(H:5, M:3, Trunc(S + 0.5):3, R:11:6);
end;
(* ------------------------------------------------------------------------ *)

begin (* Comet *)

  GetElm(T0, Q, ECC, PQR, TEQX0); (* read orbital elements *)
  GetEph(T1, DT, T2, TEQX); (* read period of time and equinox *)

  NLINE := 0;

  PrecMatEcl(TEQX0, TEQX, A); (* calculate precession matrix *)

  T := T1;

  repeat

    (* date *)

    MODJD := T * 36525.0 + 51544.5;
    CalDat(MODJD, DAY, MONTH, YEAR, HOUR);

    (* ecliptic coordinates of the sun, equinox TEQX *)

    SunPos(T, LS, BS, RS);
    Cart(RS, BS, LS, XS, YS, ZS);
    PrecMatEcl(T, TEQX, ASI);
    PrecArt(ASI, XS, YS, ZS);

    (* heliocentric ecliptic coordinates of the comet *)

    Kepler(T0, T, Q, ECC, PQR, X, Y, Z, VX, VY, VZ);
    PrecArt(A, X, Y, Z);
    PrecArt(A, VX, VY, VZ);
    Polar(X, Y, Z, R, B, L);

    (* geometric geocentric coordinates of the comet *)

    X := X + XS;
    Y := Y + YS;
    Z := Z + ZS;
    DELTA0 := SQRT(X * X + Y * Y + Z * Z);

    (* first order correction for light travel time *)

    FAC := 0.00578 * DELTA0;
    X := X - FAC * VX;
    Y := Y - FAC * VY;
    Z := Z - FAC * VZ;
    Ecl2Equ(TEQX, X, Y, Z);
    Polar(X, Y, Z, DELTA, Dec, Ra);
    Ra := Ra / 15.0;

    (* output *)

    write(YEAR:4, '/', MONTH:2, '/', DAY:2, HOUR:6:1);
    write(LS:7:1, L:7:1, B:6:1, R:7:3);
    WrtLBR(Ra, Dec, DELTA0);
    NLINE := NLINE + 1;
    if (NLINE MOD 5) = 0 then
      writeln;

    (* next time step *)

    T := T + DT;

  until (T2 < T);

end.
