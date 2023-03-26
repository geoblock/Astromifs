(* --------------------------------------------------------------------------- *)
(* planrise *)
(* planetary and solar rising and setting times *)
(* --------------------------------------------------------------------------- *)

program planrise(Input, Output);

{$APPTYPE CONSOLE}

uses
  Astro.Matlib,
  Astro.Timlib,
  Astro.Sphlib,
  Astro.Pnulib,
  Astro.Planets;

const
  J2000 = 0.0; (* Standard epoch J2000 *)
  SID = 0.9972696; (* Conversion sidereal/solar time *)
  SIN_H0P = -9.890038E-3; (* sin(-34'); altitude value for planets *)
  SIN_H0S = -1.45439E-2; (* sin(-50'); altitude value for the sun *)

type
  EVENT_TYPE = (RISING, TRANSIT, SETTING);
  STATE_TYPE = (NEVER_RISES, OK, CIRCUMPOLAR);

var
  Planet: PlanetType;
  EVENT: EVENT_TYPE;
  STATE: STATE_TYPE;
  MJD0, ZT, ZT0, D_ZT, LST, D_TAU: Double;
  LST_0H, SDA, LAMBDA, PHI, ZONE: Double;
  T, DEC, RA, SIN_H0: Double;
  SIN_PHI, COS_PHI: Double;
  RA_0H, RA_24H, DEC_0H, DEC_24H: Double;
  CNT: integer;
  PMAT: Double33;

  (* --------------------------------------------------------------------------- *)
  (* *)
  (* Formatted output *)
  (* *)
  (* --------------------------------------------------------------------------- *)

procedure WHM(UT: Double);
var
  H, M: integer;
begin
  UT := Trunc(UT * 60.0 + 0.5) / 60.0; (* rounding to 1 min *)
  H := Trunc(UT);
  M := Trunc(60.0 * (UT - H) + 0.5);
  if H < 10 then
  begin
    write(' 0');
    write(H:1, ':')
  end
  else
    write(H:3, ':');
  if M < 10 then
  begin
    write('0');
    write(M:1, ' ')
  end
  else
    write(M:2, ' ')
end;

(* --------------------------------------------------------------------------- *)
(* *)
(* Reduce (time) argument X to the interval [0..24] *)
(* *)
(* --------------------------------------------------------------------------- *)

function MOD24(X: Double): Double;
begin
  if X >= 0.0 then
    MOD24 := X - 24.0 * Trunc(X / 24.0)
  else
    MOD24 := X - 24.0 * Trunc(X / 24.0) + 24.0;
end;

(* --------------------------------------------------------------------------- *)
(* *)
(* Read time and geographical coordinates *)
(* *)
(* --------------------------------------------------------------------------- *)

procedure GetInput(var DATE, LAMBDA, PHI, ZONE: Double);

var
  D, M, Y: integer;

begin

  writeln;
  writeln(' planrise: planetary and solar rising and setting times ');
  writeln('                   version 93/07/01                     ');
  writeln('      (c) 1993 Thomas Pfleger, Oliver Montenbruck       ');
  writeln;
  write(' Date (yyyy mm dd)                    ... ');
  readln(Y, M, D);
  writeln;
  write(' Observing site: longitude (l<0 east) ... ');
  readln(LAMBDA);
  write('                 latitude             ... ');
  readln(PHI);
  write('                 local time - UT (h)  ... ');
  readln(ZONE);
  writeln;

  DATE := MJD(D, M, Y, -ZONE) (* MJD for 0h local time *)

end; (* GetInput *)

(* --------------------------------------------------------------------------- *)
(* *)
(* PLAN_RA_DEC: geocentric equatorial planetary coordinates *)
(* *)
(* Planet: Planet for which the coordinates are computed. A call with *)
(* Planet = Earth yields the geocentric coordinartes of the sun. *)
(* T     : time in julian centuries since J2000 *)
(* RA    : Right ascension in deg [0..360] *)
(* DEC   : Declination in deg *)
(* *)
(* Note: *)
(* This procedure uses the globally defined precesion matrix PMAT *)
(* for transformation from J2000 to the required equinox. *)
(* *)
(* --------------------------------------------------------------------------- *)

procedure PLAN_RA_DEC(Planet: PlanetType; T: Double; var RA, DEC: Double);
var
  XP, YP, ZP, XE, YE, ZE, R: Double;
begin
  if (Planet <> Earth) then
  begin

    POSITION(Earth, T, XE, YE, ZE);
    ECLEQU(J2000, XE, YE, ZE);

    (* Determine geocentric geometric planetary position *)
    POSITION(Planet, T, XP, YP, ZP);
    ECLEQU(J2000, XP, YP, ZP);

    XP := XP - XE;
    YP := YP - YE;
    ZP := ZP - ZE;

    PRECART(PMAT, XP, YP, ZP);

    (* Right ascension, declination and distance of the planet *)
    Polar(XP, YP, ZP, R, DEC, RA);
  end
  else

  (* Compute geocentric equatorial coordinates of the sun *)
  begin
    POSITION(Earth, T, XE, YE, ZE);
    ECLEQU(J2000, XE, YE, ZE);
    PRECART(PMAT, XE, YE, ZE);

    Polar(-XE, -YE, -ZE, R, DEC, RA);
  end;
end; (* PLAN_RA_DEC *)

(* --------------------------------------------------------------------------- *)

begin (* planrise main program *)

  (* Read input data *)
  GetInput(MJD0, LAMBDA, PHI, ZONE);

  SIN_PHI := SN(PHI);
  COS_PHI := CS(PHI);

  (* Compute local sidereal time at 0h local time *)
  LST_0H := LMST(MJD0, LAMBDA);

  (* Precession matrix (J2000 -> mean equator and equinox of date) *)
  (* for equatorial coordinates *)
  T := (MJD0 - 51544.5) / 36525.0;
  PMATEQU(J2000, T, PMAT);

  (* Compute and print rising and setting times *)
  writeln('rise':18, 'culmination':16, 'set':9);
  writeln;

  for Planet := Mercury to Pluto do
  begin

    (* Compute geocentr. planetary position at 0h and 24h local time *)
    T := (MJD0 - 51544.5) / 36525.0;
    PLAN_RA_DEC(Planet, T, RA_0H, DEC_0H);
    T := (MJD0 + 1.0 - 51544.5) / 36525.0;
    PLAN_RA_DEC(Planet, T, RA_24H, DEC_24H);

    (* Generate continuous right ascension values in case of jumps *)
    (* between 0h and 24h *)
    if (RA_0H - RA_24H) > 180.0 then
      RA_24H := RA_24H + 360.0;
    if (RA_0H - RA_24H) < -180.0 then
      RA_0H := RA_0H + 360.0;

    case Planet of
      Mercury:
        write(' Mercury ');
      Venus:
        write(' Venus   ');
      Earth:
        write(' Sun     ');
      Mars:
        write(' Mars    ');
      Jupiter:
        write(' Jupiter ');
      Saturn:
        write(' Saturn  ');
      Uranus:
        write(' Uranus  ');
      Neptune:
        write(' Neptune ');
      Pluto:
        write(' Pluto   ')
    end;
    write(' ':3);

    EVENT := RISING;
    STATE := OK;

    while ((EVENT <= SETTING) and (STATE = OK)) do
    begin
      ZT0 := 12.0; (* Starting value 12h local time *)
      CNT := 0;

      repeat
        (* Linear interpolation of planetary position *)
        RA := RA_0H + (ZT0 / 24.0) * (RA_24H - RA_0H);
        DEC := DEC_0H + (ZT0 / 24.0) * (DEC_24H - DEC_0H);

        (* Compute semi-diurnal arc (in deg) *)
        if Planet <> Earth then
          SIN_H0 := SIN_H0P
        else
          SIN_H0 := SIN_H0S;

        SDA := (SIN_H0 - SN(DEC) * SIN_PHI) / (CS(DEC) * COS_PHI);

        if (abs(SDA) < 1.0) then
        begin
          SDA := ACS(SDA);
          STATE := OK
        end
        else (* Test for circumpolar motion or invisibility *)
          if (PHI >= 0.0) then
            if DEC > (90.0 - PHI) then
              STATE := CIRCUMPOLAR
            else
              STATE := NEVER_RISES
          else if DEC < (-90.0 - PHI) then
            STATE := CIRCUMPOLAR
          else
            STATE := NEVER_RISES;

        (* Improved times for rising, culmination and setting *)
        if STATE = OK then
        begin
          LST := LST_0H + ZT0 / SID; (* Sidereal time at univ. time ZT0 *)
          case EVENT of
            RISING:
              D_TAU := (LST - RA / 15.0) + SDA / 15.0;
            TRANSIT:
              D_TAU := (LST - RA / 15.0);
            SETTING:
              D_TAU := (LST - RA / 15.0) - SDA / 15.0;
          end;
          D_ZT := SID * (MOD24(D_TAU + 12.0) - 12.0);
          ZT := ZT0 - D_ZT;
          ZT0 := ZT;
          CNT := CNT + 1
        end;

      until ((abs(D_ZT) <= 0.008) or (CNT > 10) or (STATE <> OK));

      (* Print result *)
      if STATE = OK then
      begin
        WHM(ZT);
        write(' ':6);
      end
      else
        case STATE of
          NEVER_RISES:
            write('-------- always invisible -------');
          CIRCUMPOLAR:
            write('--------- always visible ---------')
        end;

      EVENT := SUCC(EVENT);

    end;
    writeln;
  end; (* for Planet... *)

  writeln;
  write(' all times in local time ( = UT ');
  if ZONE >= 0 then
    write('+');
  writeln(ZONE:4:1, 'h )');

end.
