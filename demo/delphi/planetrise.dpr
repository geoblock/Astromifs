program planetrise(Input, Output);

(* --------------------------------------------------------------------------- *)
(* Planetary and solar rising and setting times *)
(* --------------------------------------------------------------------------- *)

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.Time,
  Apc.Spheric,
  Apc.PrecNut,
  Apc.Planets;

const
  J2000 = 0.0; // Standard epoch J2000
  SID = 0.9972696; // Conversion sidereal/solar time
  SIN_H0P = -9.890038E-3; // sin(-34'); altitude value for planets
  SIN_H0S = -1.45439E-2; // sin(-50'); altitude value for the sun

type
  EventType = (etRising, etTransit, etSetting);
  StateType = (stNeverRises, stOK, stCircumPolar);

var
  Planet: PlanetType;
  EVENT: EventType;
  STATE: StateType;
  MJD0, ZT, ZT0, D_ZT, LST, D_TAU: Double;
  LST_0H, SDA, LAMBDA, PHI, ZONE: Double;
  T, Dec, Ra, SIN_H0: Double;
  SIN_PHI, COS_PHI: Double;
  RA_0H, RA_24H, DEC_0H, DEC_24H: Double;
  CNT: integer;
  PMAT: Double33;

  (* --------------------------------------------------------------------------- *)
  (* Formatted output *)
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
(* Reduce (time) argument X to the interval [0..24] *)
(* --------------------------------------------------------------------------- *)

function Mod24(X: Double): Double;
begin
  if X >= 0.0 then
    Mod24 := X - 24.0 * Trunc(X / 24.0)
  else
    Mod24 := X - 24.0 * Trunc(X / 24.0) + 24.0;
end;

(* --------------------------------------------------------------------------- *)
(* Read time and geographical coordinates *)
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
(* Ra    : Right ascension in deg [0..360] *)
(* Dec   : Declination in deg *)
(* *)
(* Note: *)
(* This procedure uses the globally defined precesion matrix PMAT *)
(* for transformation from J2000 to the required equinox. *)
(* *)
(* --------------------------------------------------------------------------- *)

procedure PLAN_RA_DEC(Planet: PlanetType; T: Double; var Ra, Dec: Double);
var
  XP, YP, ZP, XE, YE, ZE, R: Double;
begin
  if (Planet <> Earth) then
  begin

    PlanetPos(Earth, T, XE, YE, ZE);
    Ecl2Equ(J2000, XE, YE, ZE);

    (* Determine geocentric geometric planetary position *)
    PlanetPos(Planet, T, XP, YP, ZP);
    Ecl2Equ(J2000, XP, YP, ZP);

    XP := XP - XE;
    YP := YP - YE;
    ZP := ZP - ZE;

    PrecArt(PMAT, XP, YP, ZP);

    (* Right ascension, declination and distance of the planet *)
    Polar(XP, YP, ZP, R, Dec, Ra);
  end
  else

  // Compute geocentric equatorial coordinates of the sun
  begin
    PlanetPos(Earth, T, XE, YE, ZE);
    Ecl2Equ(J2000, XE, YE, ZE);
    PrecArt(PMAT, XE, YE, ZE);

    Polar(-XE, -YE, -ZE, R, Dec, Ra);
  end;
end; // PLAN_RA_DEC

//---------------------------------------------------------------------------

begin // planetrise main program

  // Read input data
  GetInput(MJD0, LAMBDA, PHI, ZONE);

  SIN_PHI := SN(PHI);
  COS_PHI := CS(PHI);

  // Compute local sidereal time at 0h local time
  LST_0H := LMST(MJD0, LAMBDA);

  (* Precession matrix (J2000 -> mean equator and equinox of date)
    for equatorial coordinates *)
  T := (MJD0 - 51544.5) / 36525.0;
  PrecMatEqu(J2000, T, PMAT);

  // Compute and print rising and setting times
  writeln('rise':18, 'culmination':16, 'set':9);
  writeln;

  for Planet := Mercury to Pluto do
  begin
    // Compute geocentr. planetary position at 0h and 24h local time
    T := (MJD0 - 51544.5) / 36525.0;
    PLAN_RA_DEC(Planet, T, RA_0H, DEC_0H);
    T := (MJD0 + 1.0 - 51544.5) / 36525.0;
    PLAN_RA_DEC(Planet, T, RA_24H, DEC_24H);

    (* Generate continuous right ascension values in case of jumps
      between 0h and 24h *)
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

    EVENT := etRising;
    STATE := stOK;

    while ((EVENT <= etSetting) and (STATE = stOK)) do
    begin
      ZT0 := 12.0; // Starting value 12h local time
      CNT := 0;

      repeat
        // Linear interpolation of planetary position
        Ra := RA_0H + (ZT0 / 24.0) * (RA_24H - RA_0H);
        Dec := DEC_0H + (ZT0 / 24.0) * (DEC_24H - DEC_0H);

        // Compute semi-diurnal arc (in deg)
        if Planet <> Earth then
          SIN_H0 := SIN_H0P
        else
          SIN_H0 := SIN_H0S;

        SDA := (SIN_H0 - SN(Dec) * SIN_PHI) / (CS(Dec) * COS_PHI);

        if (abs(SDA) < 1.0) then
        begin
          SDA := ACS(SDA);
          STATE := stOK
        end
        else // Test for circumpolar motion or invisibility
          if (PHI >= 0.0) then
            if Dec > (90.0 - PHI) then
              STATE := stCircumPolar
            else
              STATE := stNeverRises
          else if Dec < (-90.0 - PHI) then
            STATE := stCircumPolar
          else
            STATE := stNeverRises;

        // Improved times for rising, culmination and setting
        if STATE = stOK then
        begin
          LST := LST_0H + ZT0 / SID; // Sidereal time at univ. time ZT0
          case EVENT of
            etRising:
              D_TAU := (LST - Ra / 15.0) + SDA / 15.0;
            etTransit:
              D_TAU := (LST - Ra / 15.0);
            etSetting:
              D_TAU := (LST - Ra / 15.0) - SDA / 15.0;
          end;
          D_ZT := SID * (Mod24(D_TAU + 12.0) - 12.0);
          ZT := ZT0 - D_ZT;
          ZT0 := ZT;
          CNT := CNT + 1
        end;

      until ((abs(D_ZT) <= 0.008) or (CNT > 10) or (STATE <> stOK));

      // Print result
      if STATE = stOK then
      begin
        WHM(ZT);
        write(' ':6);
      end
      else
        case STATE of
          stNeverRises:
            write('-------- always invisible -------');
          stCircumPolar:
            write('--------- always visible ---------')
        end;

      EVENT := SUCC(EVENT);

    end;
    writeln;
  end; // for Planet...

  writeln;
  write(' all times in local time ( = UT ');
  if ZONE >= 0 then
    write('+');
  writeln(ZONE:4:1, 'h )');

end.
