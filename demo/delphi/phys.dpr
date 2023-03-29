program Phys(Input, Output);

(* ----------------------------------------------------------------------- *)
(* Physical ephemeris of the major planets and the Sun *)
(* ----------------------------------------------------------------------- *)

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.Time,
  Apc.Spheric,
  Apc.PrecNut,
  Apc.Planets,
  Apc.Physic;

const
  AU = 149597870.0; (* 1 AU in km *)
  C_LIGHT = 173.14; (* Velocity of light ( AU/d ) *)
  J2000 = 0.0; (* Reference epoch J2000 *)

var
  Planet: PlanetType;
  DAY, MONTH, YEAR: integer;
  HOUR, T, T0: Double;
  X, Y, Z, XE, YE, ZE: Double;
  XX, YY, ZZ, R, DELTA: Double;
  R_EQU, F, D_EQU: Double;
  DD, DM: integer;
  DS: Double;
  L1, L2, L3, B, BSUN, D: Double;
  A0, D0, W1, W2, W3: Double;
  SENSE: RotationType;
  AX, AY, AZ, POSAX: Double;
  ELONG, PHI, K, MAG: Double;
  LSUN, DSUN, POSSUN: Double;
  PMAT: Double33;

begin
  (* Print header and read desired date *)

  writeln;
  writeln(' Phys: physical ephemeris of the planets and the Sun');
  writeln('                 Version 93/07/01                     ');
  writeln('     (c) 1993 Thomas Pfleger, Oliver Montenbruck      ');
  writeln;
  write(' Date (yyyy mm dd hh.hhh)  ');
  Readln(YEAR, MONTH, DAY, HOUR);
  writeln;

  T := (MJD(DAY, MONTH, YEAR, HOUR) - 51544.5) / 36525.0;

  write('D':13, 'V':7, 'i':7, 'PA(S)':9, 'PA(A)':8);
  writeln('L(I)':8, 'L(II)':7, 'L(III)':7, 'B':6);
  writeln('"':13, 'mag':8, 'o':6, 'o':7, 'o':8, 'o':8, 'o':7, 'o':7, 'o':8);

  (* Precession matrix (J2000 -> mean equinox of date) *)
  (* for equatorial coordinates *)

  PrecMatEqu(J2000, T, PMAT);

  (* Equatorial coordinates of the Earth, mean equinox of J2000 *)
  PlanetPos(Earth, T, XE, YE, ZE);
  Ecl2Equ(J2000, XE, YE, ZE);

  (* Physical ephemeris of the planets *)
  for Planet := Mercury to Pluto do
    if Planet <> Earth then
    begin
      (* Compute the planet's geocentric geometric position and the *)
      (* light time (in days) *)
      PlanetPos(Planet, T, X, Y, Z);
      Ecl2Equ(J2000, X, Y, Z);
      DELTA := SQRT((X - XE) * (X - XE) + (Y - YE) * (Y - YE) + (Z - ZE) * (Z - ZE));
      T0 := T - (DELTA / C_LIGHT) / 36525.0;

      (* Compute the antedated planetary position at emission of light *)
      (* (i.e. apply a first order light time correction) *)

      PlanetPos(Planet, T0, X, Y, Z);
      Ecl2Equ(J2000, X, Y, Z);

      (* Light time corrected geocentric coordinates *)
      XX := X - XE;
      YY := Y - YE;
      ZZ := Z - ZE;

      (* Compute apparent equatorial diameter (in ") *)
      Shape(Planet, R_EQU, F);
      D_EQU := 3600.0 * 2.0 * ASN(R_EQU / (DELTA * AU));

      (* Compute right ascension and declination of the axis of *)
      (* rotation with respect to the equator and equinox of J2000; *)
      (* compute orientation of the prime meridian *)

      Orient(Planet, SYS_I, T0, A0, D0, W1, SENSE);
      Orient(Planet, SYS_II, T0, A0, D0, W2, SENSE);
      Orient(Planet, SYS_III, T0, A0, D0, W3, SENSE);

      (* Compute planetocentric longitude and latitude of the Earth *)
      Rotation(XX, YY, ZZ, A0, D0, W1, SENSE, F, AX, AY, AZ, L1, B, D);
      Rotation(XX, YY, ZZ, A0, D0, W2, SENSE, F, AX, AY, AZ, L2, B, D);
      Rotation(XX, YY, ZZ, A0, D0, W3, SENSE, F, AX, AY, AZ, L3, B, D);

      (* Compute planetocentric longitude and latitude of the Sun *)
      Rotation(X, Y, Z, A0, D0, W1, SENSE, F, AX, AY, AZ, LSUN, BSUN, DSUN);

      (* Compute illumination and apparent magnitude *)
      Illum(X, Y, Z, XE, YE, ZE, R, DELTA, ELONG, PHI, K);
      MAG := Bright(Planet, R, DELTA, PHI, DSUN, LSUN - L1);

      (* Compute position angles of the axis of rotation and of the *)
      (* Sun with respect to the mean equinox of date *)
      PrecArt(PMAT, X, Y, Z);
      PrecArt(PMAT, XX, YY, ZZ);
      PrecArt(PMAT, AX, AY, AZ);

      POSAX := PosAng(XX, YY, ZZ, AX, AY, AZ);
      POSSUN := PosAng(XX, YY, ZZ, -X, -Y, -Z);

      (* Print results *)
      case Planet of
        Mercury:
          write(' Mercury  ');
        Venus:
          write(' Venus    ');
        Mars:
          write(' Mars     ');
        Jupiter:
          write(' Jupiter  ');
        Saturn:
          write(' Saturn   ');
        Uranus:
          write(' Uranus   ');
        Neptune:
          write(' Neptune  ');
        Pluto:
          write(' Pluto    ')
      end;

      write(D_EQU:5:2, MAG:6:1, PHI:7:1);
      write(POSSUN:8:2, POSAX:8:2);

      case Planet of
        Mercury, Venus, Mars, Pluto:
          write(L1:8:2, ' ':14);
        Jupiter:
          write(L1:8:2, L2:7:2, L3:7:2);
        Saturn:
          write(L1:8:2, L3:14:2);
        Uranus, Neptune:
          write(L3:22:2)
      end;

      writeln(B:8:2);
    end;
  (* Physical ephemeris of the Sun *)

  (* Compute light time corrected equatorial coordinates of *)
  (* the Earth with respect to the equator and equinox of J2000 *)

  DELTA := SQRT(XE * XE + YE * YE + ZE * ZE);
  T0 := T - (DELTA / C_LIGHT) / 36525.0;
  PlanetPos(Earth, T0, XE, YE, ZE);
  Ecl2Equ(J2000, XE, YE, ZE);

  (* Right ascension and declination of the Sun's axis (J2000), *)
  (* orientation of the prime meridian and equatorial radius (km) *)

  A0 := 285.96;
  D0 := 63.96;
  W1 := 84.11 + 14.1844000 * 36525.0 * T;
  W1 := W1 / 360.0;
  W1 := 360.0 * (W1 - Trunc(W1));
  R_EQU := 696000.0;

  (* Compute heliographic coordinates of the Earth *)
  Rotation(-XE, -YE, -ZE, A0, D0, W1, Retrograde, 0.0, AX, AY, AZ, L1, B, B);

  (* Compute position angle of the axis of rotation *)
  (* with respect to the mean equinox of date *)
  PrecArt(PMAT, XE, YE, ZE);
  PrecArt(PMAT, AX, AY, AZ);
  POSAX := PosAng(-XE, -YE, -ZE, AX, AY, AZ);

  (* Express position angle of the Sun's axis within -180..180 deg *)
  if POSAX > 180.0 then
    POSAX := POSAX - 360.0;

  (* Compute apparent equatorial diameter (in ") *)
  D_EQU := 2.0 * ASN(R_EQU / (DELTA * AU));

  // Output results
  DMS(D_EQU, DD, DM, DS);
  writeln('''':10, '"':3, 'o':29, 'o':8, 'o':22);
  writeln(' Sun  ', DM:3, DS:6:2, POSAX:29:2, L1:8:2, B:22:2)

end.
