(* ----------------------------------------------------------------------- *)
(* PlanetPos *)
(* heliocentric and geocentric planetary positions *)
(* ----------------------------------------------------------------------- *)

program PlanetPos(Input, Output);

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric,
  Apc.Sun,
  Apc.Planets,
  Apc.Time;

const
  J2000 = 0.0;
  B1950 = -0.500002108;

var
  DAY, MONTH, YEAR, IPLAN, IMODE, K: integer;
  HOUR, MODJD, T, TEQX: Double;
  X, Y, Z, XP, YP, ZP, XS, YS, ZS: Double;
  L, B, R, LS, BS, RS, Ra, Dec, DELTA, DELTA0: Double;
  A: Double33;
  CharMode: Char;

  (* ----------------------------------------------------------------------- *)
  (* PRINTOUT: print coordinates of one planet *)
  (* ----------------------------------------------------------------------- *)
procedure PRINTOUT(IPLAN: integer; L, B, R, Ra, Dec, DELTA: Double);
var
  H, M: integer;
  S: Double;
begin
  DMS(L, H, M, S);
  write(H:3, M:3, S:5:1);
  DMS(B, H, M, S);
  write(H:4, M:3, S:5:1);
  if IPLAN < 4 then
    write(R:11:6)
  else
    write(R:10:5, ' ');
  DMS(Ra, H, M, S);
  write(H:4, M:3, S:6:2);
  DMS(Dec, H, M, S);
  write(H:4, M:3, S:5:1);
  if IPLAN < 4 then
    write(DELTA:11:6)
  else
    write(DELTA:10:5, ' ');
end;
(* ----------------------------------------------------------------------- *)

begin (* PlanetPos *)

  writeln;
  writeln('   PlanetPos: geocentric and heliocentric planetary positions ');
  writeln('                       version 93/07/01                     ');
  writeln('         (c) 1993 Thomas Pfleger, Oliver Montenbruck        ');
  writeln;

  repeat

    writeln;
    writeln(' (J) J2000 astrometric       (B) B1950 astrometric  ');
    writeln(' (A) apparent coordinates    (E) end                ');
    writeln;
    write(' enter option: ');
    Readln(CharMode);
    writeln;

    if CharMode in ['A', 'a', 'J', 'j', 'B', 'b'] then

    begin

      (* read date *)

      write(' date (year month day hour) ?     ');
      Readln(YEAR, MONTH, DAY, HOUR);
      writeln;
      writeln;
      writeln;
      writeln;
      MODJD := MJD(DAY, MONTH, YEAR, HOUR);
      T := (MODJD - 51544.5) / 36525.0;
      write(' date:  ', YEAR:4, '/', MONTH:2, '/', DAY:2, ' ', HOUR:5:1, '(ET)');
      write('JD:':6, (MODJD + 2400000.5):12:3, 'equinox ':18);
      case CharMode of
        'A', 'a':
          writeln('of date');
        'J', 'j':
          writeln('J2000');
        'B', 'b':
          writeln('B1950');
      end;
      writeln;

      (* header *)

      write(' ':10, 'l':6, 'b':12, 'r':11);
      writeln(' ':7, 'Ra':5, 'Dec':13, 'delta':13);
      write(' ':9, '   o  ''  "', ' ':3, '  o  ''  "', ' ':6, 'AU', ' ':4);
      writeln(' ':2, '  h  m  s', ' ':4, '  o  ''  "', ' ':6, 'AU');

      (* ecliptic coordinates of the sun, Equinox T *)

      SunPos(T, LS, BS, RS);

      (* planetary coordinates; include Pluto between 1890 and 2100 *)

      if ((-1.1 < T) and (T < +1.0)) then
        K := 9
      else
        K := 8;

      for IPLAN := 0 to K do

      begin

        (* heliocentric ecliptic coordinates of the planet *)

        case IPLAN of
          1:
            MercuryPos(T, L, B, R);
          2:
            VenusPos(T, L, B, R);
          4:
            MarsPos(T, L, B, R);
          5:
            JupiterPos(T, L, B, R);
          6:
            SaturnPos(T, L, B, R);
          7:
            UranusPos(T, L, B, R);
          8:
            NeptunePos(T, L, B, R);
          9:
            PlutoPos(T, L, B, R);
          0:
            begin
              L := 0.0;
              B := 0.0;
              R := 0.0;
            end;
          3:
            begin
              L := LS + 180.0;
              B := -BS;
              R := RS;
            end;
        end;

        (* geocentric ecliptic coordinates (light-time corrected) *)

        if CharMode in ['A', 'a'] then
          IMODE := 2
        else
          IMODE := 1;
        GeoCentric(T, L, B, R, LS, BS, RS, IPLAN, IMODE, XP, YP, ZP, XS, YS, ZS, X, Y, Z, DELTA0);

        (* precession, equatorial coordinates, nutation *)

        case CharMode of
          'J', 'j':
            TEQX := J2000;
          'B', 'b':
            TEQX := B1950;
        end;

        if CharMode in ['A', 'a'] then
        begin
          Ecl2Equ(T, X, Y, Z);
          NutEqu(T, X, Y, Z);
        end
        else
        begin
          PrecMatEcl(T, TEQX, A);
          PrecArt(A, XP, YP, ZP);
          PrecArt(A, X, Y, Z);
          Ecl2Equ(TEQX, X, Y, Z);
        end;

        (* spherical coordinates *)

        Polar(XP, YP, ZP, R, B, L);
        Polar(X, Y, Z, DELTA, Dec, Ra);
        Ra := Ra / 15.0;

        (* output *)

        case IPLAN of
          0:
            write(' Sun     ');
          1:
            write(' Mercury ');
          2:
            write(' Venus   ');
          3:
            write(' Earth   ');
          4:
            write(' Mars    ');
          5:
            write(' Jupiter ');
          6:
            write(' Saturn  ');
          7:
            write(' Uranus  ');
          8:
            write(' Neptune ');
          9:
            write(' Pluto   ');
        end;

        PRINTOUT(IPLAN, L, B, R, Ra, Dec, DELTA0);
        writeln;

      end;

      writeln;
      writeln(' l,b,r:   heliocentric ecliptic (geometric) ');
      write(' Ra,Dec:  geocentric equatorial ');
      if CharMode in ['A', 'a'] then
        writeln('(apparent)')
      else
        writeln('(astrometric)');
      writeln(' delta:   geocentric distance   (geometric)');
      writeln;

    end;

  until CharMode in ['E', 'e']

  end. (* PlanetPos *)
