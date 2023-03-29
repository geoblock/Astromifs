program OrbDet(Input, Output, OrbInput, OrbOutput);

(* ----------------------------------------------------------------------- *)
(* Gaussian orbit determination from three observations *)
(* using the abbreviated method of Bucerius *)
(* ----------------------------------------------------------------------- *)

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric,
  Apc.Sun,
  Apc.Time,
  Apc.Kepler;

type
  CHAR80 = array [1 .. 80] of Char;

var
  TEQX: Double;
  TP, Q, ECC, INC, LAN, AOP: Double;
  JD0: Double3;
  RSUN, E: Mat3X;
  HEADER: CHAR80;
  OrbInput, OrbOutput: TEXT;

  (* -----------------------------------------------------------------------
   START: reads the input file and preprocesses the observational data output:
   RSUN:  matrix of three Sun position vectors in ecliptic coordinates
   E:     matrix of three observation direction unit vectors
   JD:    julian date of the three observation times
   TEQX:  equinox of RSUN and E (in Julian centuries since J2000)
   ----------------------------------------------------------------------- *)

procedure START(var HEADER: CHAR80; var RSUN, E: Mat3X; var JD0: Double3; var TEQX: Double);
var
  DAY, MONTH, YEAR, D, M, I: Integer;
  UT, S, DUMMY: Double;
  EQX0, EQX, TEQX0: Double;
  LS, BS, RS, LP, BP, Ra, Dec, T: Double3;
  A, ASI: Double33;
  OrbInput: TEXT;

begin
  // open input file
  Assign(OrbInput, 'ORBINP.DAT');
  Reset(OrbInput);

  // read data from file OrbInput
  for I := 1 to 80 do // header
    if not(EOLN(OrbInput)) then
      read(OrbInput, HEADER[I])
    else
      HEADER[I] := ' ';
  readln(OrbInput);
  for I := 1 to 3 do // 3 observations
  begin
    Read(OrbInput, YEAR, MONTH, DAY, UT); (* date *)
    Read(OrbInput, D, M, S);
    Ddd(D, M, S, Ra[I]); (* Ra *)
    readln(OrbInput, D, M, S);
    Ddd(D, M, S, Dec[I]); (* Dec *)
    Ra[I] := Ra[I] * 15.0;
    JD0[I] := 2400000.5 + MJD(DAY, MONTH, YEAR, UT);
    T[I] := (JD0[I] - 2451545.0) / 36525.0;
  end;
  writeln;
  readln(OrbInput, EQX0);
  TEQX0 := (EQX0 - 2000.0) / 100.0; // equinox

  // desired equinox of the orbital elements

  read(OrbInput, EQX);
  TEQX := (EQX - 2000.0) / 100.0;

  // calculate initial data of the orbit determination

  PrecMatEcl(TEQX0, TEQX, A);
  for I := 1 to 3 do
  begin
    Cart(1.0, Dec[I], Ra[I], E[I, X], E[I, Y], E[I, Z]);
    Equ2Ecl(TEQX0, E[I, X], E[I, Y], E[I, Z]);
    PrecArt(A, E[I, X], E[I, Y], E[I, Z]);
    Polar(E[I, X], E[I, Y], E[I, Z], DUMMY, BP[I], LP[I]);
    PrecMatEcl(T[I], TEQX, ASI);
    SunPos(T[I], LS[I], BS[I], RS[I]);
    Cart(RS[I], BS[I], LS[I], RSUN[I, X], RSUN[I, Y], RSUN[I, Z]);
    PrecArt(ASI, RSUN[I, X], RSUN[I, Y], RSUN[I, Z]);
  end;

  writeln('   ORBDET: orbit determination from three observations ');
  writeln('                     version 93/07/01                  ');
  writeln('        (c) 1993 Thomas Pfleger, Oliver Montenbruck    ');
  writeln;
  writeln;
  writeln(' Summary of orbit determination ');
  writeln;
  write('  ');
  for I := 1 to 78 do
    write(HEADER[I]);
  writeln;
  writeln;
  writeln(' Initial data (ecliptic geocentric coordinates (in deg))');
  writeln;
  writeln('  Julian Date        ', JD0[1]:12:2, JD0[2]:12:2, JD0[3]:12:2);
  writeln('  Solar longitude    ', LS[1]:12:2, LS[2]:12:2, LS[3]:12:2);
  writeln('  Planet/Comet Longitude', LP[1]:9:2, LP[2]:12:2, LP[3]:12:2);
  writeln('  Planet/Comet Latitude ', BP[1]:9:2, BP[2]:12:2, BP[3]:12:2);
  writeln;
  writeln;

  Readln;
end;

(* ----------------------------------------------------------------------- *)
(* DumpElem: output of orbital elements (screen) *)
(* ----------------------------------------------------------------------- *)
procedure DumpElem(TP, Q, ECC, INC, LAN, AOP, TEQX: Double);
var
  DAY, MONTH, YEAR: Integer;
  MODJD, UT: Double;
begin
  MODJD := TP * 36525.0 + 51544.5;
  CalDat(MODJD, DAY, MONTH, YEAR, UT);
  writeln(' Orbital elements', ' (Equinox ', 'J', 100.0 * TEQX + 2000.0:8:2, ')');
  writeln;
  writeln('  Perihelion date      tp    ', YEAR:4, '/', MONTH:2, '/', DAY:2, UT:8:4, 'h', '  (JD',
    MODJD + 2400000.5:11:2, ')');
  writeln('  Perihelion distance  q[AU] ', Q:12:6);
  writeln('  Semi-major axis      a[AU] ', Q / (1 - ECC):12:6);
  writeln('  Eccentricity         e     ', ECC:12:6);
  writeln('  Inclination          i     ', INC:10:4, ' degrees');
  writeln('  Ascending node       Omega ', LAN:10:4, ' degrees');
  writeln('  Long. of perihelion  pi    ', AOP + LAN:10:4, ' degrees');
  writeln('  Arg. of perihelion   omega ', AOP:10:4, ' degrees');
  writeln;
end;

(* ----------------------------------------------------------------------- *)
(* SaveElem: output of orbital elements (file) *)
(* ----------------------------------------------------------------------- *)
procedure SaveElem(TP, Q, ECC, INC, LAN, AOP, TEQX: Double; HEADER: CHAR80);

var
  I, DAY, MONTH, YEAR: Integer;
  MODJD, UT: Double;

begin
  // open file for writing

  (* REWRITE(OrbOut); *)
  Assign(OrbOutput, 'OrbOut.DAT');
  REWRITE(OrbOutput);

  MODJD := TP * 36525.0 + 51544.5;
  CalDat(MODJD, DAY, MONTH, YEAR, UT);
  write(OrbOutput, YEAR:5, MONTH:3, (DAY + UT / 24.0):7:3, '!':6);
  writeln(OrbOutput, ' perihelion time T0 (y m d.d)  =  JD ', (MODJD + 2400000.5):12:3);
  writeln(OrbOutput, Q:12:6, '!':9, ' q  ( a =', Q / (1 - ECC):10:6, ' )');
  writeln(OrbOutput, ECC:12:6, '!':9, ' e ');
  writeln(OrbOutput, INC:10:4, '!':11, ' i ');
  writeln(OrbOutput, LAN:10:4, '!':11, ' long.asc.node ');
  writeln(OrbOutput, AOP:10:4, '!':11, ' arg.perih. ( long.per. = ', AOP + LAN:9:4, ' )');
  writeln(OrbOutput, TEQX * 100.0 + 2000.0:8:2, '!':13, ' equinox (J)');
  write(OrbOutput, '! ');
  for I := 1 to 78 do
    write(OrbOutput, HEADER[I]);

  Reset(OrbOutput); (* close file *)

end;

(* ----------------------------------------------------------------------- *)
(* RETARD: light-time correction *)
(* JD0: times of observation (t1',t2',t3') (Julian Date) *)
(* RHO: three geocentric distances (in AU) *)
(* JD:  times of light emittance (t1,t2,t3) (Julian Date) *)
(* TAU: scaled time differences *)
(* ----------------------------------------------------------------------- *)
procedure RETARD(JD0, RHO: Double3; var JD, TAU: Double3);
const
  KGAUSS = 0.01720209895;
  A = 0.00578;
var
  I: Integer;
begin
  for I := 1 to 3 do
    JD[I] := JD0[I] - A * RHO[I];
  TAU[1] := KGAUSS * (JD[3] - JD[2]);
  TAU[2] := KGAUSS * (JD[3] - JD[1]);
  TAU[3] := KGAUSS * (JD[2] - JD[1]);
end;

//----------------------------------------------------------------------
(* Gauss: iteration of the abbreviated Gauss method *)
(* RSUN: three vectors of geocentric Sun positions *)
(* E   : three unit vectors of geocentric observation directions *)
(* JD0 : three observation times (Julian Date) *)
(* TP  : time of perihelion passage (Julian centuries since J2000) *)
(* Q   : perihelion distance *)
(* ECC : eccentricity *)
(* INC : inclination *)
(* LAN : longitude of the ascending node *)
(* AOP : argument of perihelion *)
//----------------------------------------------------------------------

procedure Gauss(RSUN, E: Mat3X; JD0: Double3; var TP, Q, ECC, INC, LAN, AOP: Double);

const
  EPS_RHO = 1.0E-8;

var
  I, J: Integer;
  S: Index;
  RHOOLD, DET: Double;
  JD, RHO, N, TAU, ETA: Double3;
  DI: Vector;
  RPL: Mat3X;
  DD: Double33;

begin

  (* calculate initial approximations of n1 and n3 *)

  N[1] := (JD0[3] - JD0[2]) / (JD0[3] - JD0[1]);
  N[2] := -1.0;
  N[3] := (JD0[2] - JD0[1]) / (JD0[3] - JD0[1]);

  (* calculate matrix D and its determinant (det(D) = e3.d3) *)

  Cross(E[2], E[3], DI);
  for J := 1 to 3 do
    DD[1, J] := Dot(DI, RSUN[J]);
  Cross(E[3], E[1], DI);
  for J := 1 to 3 do
    DD[2, J] := Dot(DI, RSUN[J]);
  Cross(E[1], E[2], DI);
  for J := 1 to 3 do
    DD[3, J] := Dot(DI, RSUN[J]);
  DET := Dot(E[3], DI);

  writeln;
  writeln(' Iteration of the geocentric distances rho [AU] ');
  writeln;

  RHO[2] := 0;

  // Iterate until distance rho[2] does not change any more
  RHO[2] := 0;
  repeat
    RHOOLD := RHO[2];

    // geocentric distance rho from n1 and n3
    for I := 1 to 3 do
      RHO[I] := (N[1] * DD[I, 1] - DD[I, 2] + N[3] * DD[I, 3]) / (N[I] * DET);

    // apply light-time correction and calculate time differences
    RETARD(JD0, RHO, JD, TAU);

    // heliocentric coordinate vectors
    for I := 1 to 3 do
      for S := X to Z do
        RPL[I, S] := RHO[I] * E[I, S] - RSUN[I, S];

    // sector/triangle ratios eta[i]
    ETA[1] := Find_ETA(RPL[2], RPL[3], TAU[1]);
    ETA[2] := Find_ETA(RPL[1], RPL[3], TAU[2]);
    ETA[3] := Find_ETA(RPL[1], RPL[2], TAU[3]);

    // improvement of the sector/triangle ratios
    N[1] := (TAU[1] / ETA[1]) / (TAU[2] / ETA[2]);
    N[3] := (TAU[3] / ETA[3]) / (TAU[2] / ETA[2]);
    writeln('  rho', ' ': 16, RHO[1]: 12: 8, RHO[2]: 12: 8, RHO[3]: 12: 8);

  until (abs(RHO[2] - RHOOLD) < EPS_RHO);

  writeln;
  writeln(' Heliocentric distances [AU]:');
  writeln;
  writeln('  r  ', ' ':16, Norm(RPL[1]):12:8, Norm(RPL[2]):12:8, Norm(RPL[3]):12:8);
  writeln;
  writeln;

  // derive orbital elements from first and third observation
  Element(JD[1], JD[3], RPL[1], RPL[3], TP, Q, ECC, INC, LAN, AOP);

end;

//----------------------------------------------------------------------

begin
  START(HEADER, RSUN, E, JD0, TEQX);
  Gauss(RSUN, E, JD0, TP, Q, ECC, INC, LAN, AOP);
  DumpElem(TP, Q, ECC, INC, LAN, AOP, TEQX);
  SaveElem(TP, Q, ECC, INC, LAN, AOP, TEQX, HEADER);

  // check solution
  writeln;
  if (Dot(E[2], RSUN[2]) > 0) then
    writeln(' Warning: observation in hemisphere of conjunction;', '  possible second solution');
  if (ECC > 1.1) then
    writeln(' Warning: probably not a realistic solution (e>1.1) ');
  if ((abs(Q - 0.985) < 0.1) and (abs(ECC - 0.015) < 0.05)) then
    writeln(' Warning: probably Earth''s orbit solution');

end. // Orbdet


