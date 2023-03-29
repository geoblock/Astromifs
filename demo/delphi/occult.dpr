(* ----------------------------------------------------------------------- *)
(* Occult *)
(* prediction of stellar occultations by the Moon *)
(* ----------------------------------------------------------------------- *)

program Occult(Input, Output, OCCINP);

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric,
  Apc.Sun,
  Apc.Moon,
  Apc.Time;

const
  TOVLAP = 3.42E-6; (* 3h in julian centuries *)
  T_SEARCH = 2.737850787E-4; (* 10d in julian centuries *)
  NAME_LENGTH = 17; (* maxim.length of a star's name *)

type
  NAME_STRING = array [1 .. NAME_LENGTH] of Char;

var
  T_begin, T_END, T1, T2, TM, T_EQX, T_EPOCH: Double;
  ETDIFUT: Double;
  RA_STAR, DE_STAR: Double;
  VX, VY, VZ: Double;
  LAMBDA, PHI, RCPHI, RSPHI: Double;
  RAPOLY, DEPOLY, RPOLY: TPolynomCheb;
  PNMAT: Double33;
  OCCINP: TEXT;
  NAME: NAME_STRING;

  (* ----------------------------------------------------------------------- *)
  (* Get_Input: read desired period of time and observer's coordinates *)
  (* ----------------------------------------------------------------------- *)
procedure Get_Input(var T_begin, T_END, ETDIFUT, LAMBDA, PHI: Double);
var
  D, M, Y: integer;
  T: Double;
  VALID: Boolean;
begin
  writeln;
  writeln('          Occult: occultations of stars by the Moon       ');
  writeln('                      version 93/07/01                    ');
  writeln('        (c) 1993 Thomas Pfleger, Oliver Montenbruck       ');
  writeln;
  writeln(' Period of time for prediction of occultations ');
  write('   first date (yyyy mm dd)                   ... ');
  readln(Y, M, D);
  T_begin := (MJD(D, M, Y, 0) - 51544.5) / 36525.0;
  write('   last  date (yyyy mm dd)                   ... ');
  readln(Y, M, D);
  T_END := (MJD(D, M, Y, 0) - 51544.5) / 36525.0;
  T := (T_begin + T_END) / 2.0;
  ETminUT(T, ETDIFUT, VALID);
  if (VALID) then
    write(' Difference ET-UT (proposal:', Trunc(ETDIFUT + 0.5):3, ' sec)          ... ')
  else
    write(' Difference ET-UT (sec)                      ... ');
  readln(ETDIFUT);
  write(' Observer''s coordinates: longitude (>0 west) ... ');
  readln(LAMBDA);
  write('                         latitude            ... ');
  readln(PHI);
end;

(* ----------------------------------------------------------------------- *)
(* HEADER: print header *)
(* ----------------------------------------------------------------------- *)
procedure HEADER;
begin
  writeln;
  writeln('   Date        UT     D/R   Pos ', '    h      a     b       Star');
  writeln('              h  m  s          o', '     o     m     m           ');
end;

(* ----------------------------------------------------------------------- *)
(* GETSTAR: read star coordinates from file OCCINP and correct for *)
(* proper motion *)
(* ----------------------------------------------------------------------- *)
procedure GETSTAR(T_EPOCH, T: Double; var Ra, Dec: Double; var NAME: NAME_STRING);
var
  G, M, I: integer;
  S, PM_RA, PM_DEC: Double;
begin
  read(OCCINP, G, M, S);
  Ddd(G, M, S, Ra); (* right ascension at epoch *)
  read(OCCINP, PM_RA);
  read(OCCINP, G, M, S);
  Ddd(G, M, S, Dec); (* declination at epoch *)
  read(OCCINP, PM_DEC);
  Ra := Ra + (T - T_EPOCH) * PM_RA / 3600.0; (* proper motion right asc. *)
  Dec := Dec + (T - T_EPOCH) * PM_DEC / 3600.0; (* proper motion declination *)
  for I := 1 to NAME_LENGTH do (* name of the star *)
    if (not EOLN(OCCINP)) then
      read(OCCINP, NAME[I])
    else
      NAME[I] := ' ';
  readln(OCCINP);
  Ra := 15.0 * Ra; (* Ra in deg *)
end;

(* ----------------------------------------------------------------------- *)
(* Conjunct: *)
(* *)
(* checks whether there is a conjunction of Moon and star between TA *)
(* and TB during which the Moon's shadow hits the Earth *)
(* *)
(* TA,TB:    time interval for search of conjunction *)
(* RAPOLY,DEPOLY,RPOLY: Chebyshev coefficients for lunar coordinates *)
(* Ra,Dec:   right ascension and declination of the star (0<=Ra<=360) *)
(* CONJ:     True/False (conjunction found / no conjunction found) *)
(* T_CONJ:   time of conjunction in right ascension (0.0 if CONJ=False) *)
(* *)
(* All times are counted in Julian centuries since J2000. *)
(* The Chebyshev expansion of the lunar right ascension has to yield *)
(* values between -360 and +360 degress and cover less than one orbit. *)
(* ----------------------------------------------------------------------- *)

procedure Conjunct(TA, TB: Double; RAPOLY, DEPOLY, RPOLY: TPolynomCheb; Ra, Dec: Double;
  var CONJ: Boolean; var T_CONJ: Double);

const
  EPS = 1E-3; (* accuracy in degrees Ra *)

var
  RA_A, RA_B: Double;
  T1, T2, T_NEW, DRA1, DRA2, DRA_NEW: Double;
  DE_CON, R_CON: Double;

begin

  T_CONJ := 0.0;
  RA_A := T_Eval(RAPOLY, TA);
  RA_B := T_Eval(RAPOLY, TB);

  (* check if RA_A <= Ra <= RA_B *)
  CONJ := (RA_A <= Ra) and (Ra <= RA_B);
  if (not CONJ) then (* check again with Ra-360deg *)
  begin
    Ra := Ra - 360.0;
    CONJ := ((RA_A <= Ra) and (Ra <= RA_B));
  end;

  if CONJ then

  begin

    (* determine time of conjunction using 'regula falsi' *)
    (* ([T1,T2] always contains T_CONJ) *)

    T1 := TA;
    DRA1 := RA_A - Ra;
    T2 := TB;
    DRA2 := RA_B - Ra;
    repeat
      T_NEW := T2 - DRA2 * (T2 - T1) / (DRA2 - DRA1);
      DRA_NEW := T_Eval(RAPOLY, T_NEW) - Ra;
      if DRA_NEW > 0 then
      begin
        T2 := T_NEW;
        DRA2 := DRA_NEW
      end
      else
      begin
        T1 := T_NEW;
        DRA1 := DRA_NEW
      end;
    until (abs(DRA_NEW) < EPS);
    T_CONJ := T_NEW;

    (* check if lunar shadow hits the Earth *)

    DE_CON := T_Eval(DEPOLY, T_CONJ);
    R_CON := T_Eval(RPOLY, T_CONJ);
    CONJ := (abs(SN(DE_CON - Dec) * R_CON) < 1.5);

  end;

end;

(* ----------------------------------------------------------------------- *)
(* Shadow: *)
(* *)
(* starting from the time of conjunction the times, position angles *)
(* and longitude and latitude coefficients of disappearance and *)
(* reappearance are calculated for a specific observing site *)
(* *)
(* RAPOLY,DEPOLY,RPOLY  : Chebyshev approximations of lunar coordinates *)
(* T_CONJ_ET            : time of conjunction in right ascension *)
(* (in Julian centuries ET since J2000) *)
(* ETDIFUT              : ET-UT in sec *)
(* LAMBDA,RCPHI,RSPHI   : geocentric coordinates of the observer *)
(* RA_STAR,DE_STAR      : right ascension and declination of the star *)
(* EVENT                : True occultation takes place; False otherwise *)
(* MJD_UT_IN, MJD_UT_OUT: times of contact (Modified Julian Date UT) *)
(* POS_IN, POS_OUT      : position angles *)
(* H_IN, H_OUT          : star's altitude above the horizon *)
(* A_IN, A_OUT          : longitude coefficient *)
(* B_IN, B_OUT          : latitude coefficient *)
(* ----------------------------------------------------------------------- *)

procedure Shadow

  (RAPOLY, DEPOLY, RPOLY: TPolynomCheb; T_CONJ_ET, ETDIFUT, LAMBDA, RCPHI, RSPHI, RA_STAR,
  DE_STAR: Double; var EVENT: Boolean; var MJD_UT_IN, MJD_UT_OUT: Double;
  var POS_IN, POS_OUT, H_IN, H_OUT, A_IN, A_OUT, B_IN, B_OUT: Double);

const
  DTAB = 0.25; (* search step size in hours *)
  RANGE = 2.25; (* search interval = +/-(RANGE+DTAB) in hours *)
  K = 0.2725; (* ratio earth radius / lunar radius *)
  CENT = 876600.0; (* hours per Julian century *)
  SID = 1.0027379; (* ratio solar time / sidereal time *)

var
  I, NZ, NFOUND: integer;
  K_SQR, MJD_CONJ_UT, HOUR, F, G: Double;
  THETA_CONJ, CSDEST, SNDEST: Double;
  S_MINUS, S_0, S_PLUS, XE, YE: Double;
  Z, TIME: array [1 .. 2] of Double;

  (* FG: f-g coordinates in the fundamental plane *)
  procedure FG(HOUR: Double; var F, G: Double);
  var
    T, DEM, RM, RCDEM, RSDEM, DRAM, DRA: Double;
  begin
    T := T_CONJ_ET + HOUR / CENT;
    DEM := T_Eval(DEPOLY, T);
    RM := T_Eval(RPOLY, T);
    RCDEM := RM * CS(DEM);
    RSDEM := RM * SN(DEM);
    DRAM := T_Eval(RAPOLY, T) - RA_STAR;
    DRA := 15.0 * (THETA_CONJ + HOUR * SID) - RA_STAR;
    F := +RCDEM * SN(DRAM) - RCPHI * SN(DRA);
    G := +RSDEM * CSDEST - RCDEM * SNDEST * CS(DRAM) - RSPHI * CSDEST + RCPHI * SNDEST * CS(DRA);
  end;

(* CONTACT: position angle, altitude and longitude/latitude coeffic. *)
  procedure CONTACT(HOUR: Double; var POS, H, A, B: Double);
  var
    F, G, FF, GG, DF, DG, FAC, DRA, CDRA, SDRA: Double;
  begin
    FG(HOUR, F, G);
    FG(HOUR + DTAB, FF, GG);
    DF := (FF - F) / DTAB;
    DG := (GG - G) / DTAB;
    POS := ATN2(-F, -G);
    if POS < 0.0 then
      POS := POS + 360.0;
    FAC := 1.047 / (F * DF + G * DG);
    DRA := 15.0 * (THETA_CONJ + HOUR * SID) - RA_STAR;
    CDRA := CS(DRA);
    SDRA := SN(DRA);
    A := -FAC * RCPHI * (F * CDRA + G * SDRA * SNDEST);
    B := -FAC * (RSPHI * (F * SDRA - G * SNDEST * CDRA) - RCPHI * G * CSDEST);
    H := ASN(RSPHI * SNDEST + RCPHI * CSDEST * CDRA);
  end;

begin

  (* modified julian date and sidereal time at time of conjunction *)
  MJD_CONJ_UT := T_CONJ_ET * 36525.0 + 51544.5 - ETDIFUT / 86400.0;
  THETA_CONJ := LMST(MJD_CONJ_UT, LAMBDA);

  (* auxiliary values *)
  K_SQR := K * K;
  CSDEST := CS(DE_STAR);
  SNDEST := SN(DE_STAR);

  (* search for time of contact *)
  NFOUND := 0;
  TIME[1] := 0.0;
  TIME[2] := 0.0;
  HOUR := -RANGE - 2.0 * DTAB;
  FG(-RANGE - DTAB, F, G);
  S_PLUS := F * F + G * G - K_SQR;
  repeat
    HOUR := HOUR + 2.0 * DTAB;
    S_MINUS := S_PLUS;
    FG(HOUR, F, G);
    S_0 := F * F + G * G - K_SQR;
    FG(HOUR + DTAB, F, G);
    S_PLUS := F * F + G * G - K_SQR;
    Quad(S_MINUS, S_0, S_PLUS, XE, YE, Z[1], Z[2], NZ);
    for I := 1 to NZ do
      TIME[NFOUND + I] := HOUR + DTAB * Z[I];
    NFOUND := NFOUND + NZ;
    EVENT := (NFOUND = 2);
  until ((EVENT) or (HOUR >= RANGE));

  (* calculate details of an occultation *)
  if EVENT then
  begin
    MJD_UT_IN := MJD_CONJ_UT + TIME[1] / 24.0;
    MJD_UT_OUT := MJD_CONJ_UT + TIME[2] / 24.0;
    CONTACT(TIME[1], POS_IN, H_IN, A_IN, B_IN);
    CONTACT(TIME[2], POS_OUT, H_OUT, A_OUT, B_OUT);
  end;

end;

(* ----------------------------------------------------------------------- *)
(* DARKNESS: test for civil twilight *)
(* MODJD:     Modified Julian Date *)
(* LAMBDA:    geographic longitude (>0 west of Greenwich) *)
(* CPHI,SPHI: sine and cosine of the geographic latitude *)
(* ----------------------------------------------------------------------- *)
function DARKNESS(MODJD, LAMBDA, CPHI, SPHI: Double): Boolean;
var
  T, Ra, Dec, TAU, SIN_HSUN: Double;
begin
  T := (MODJD - 51544.5) / 36525.0;
  MiniSun(T, Ra, Dec);
  TAU := 15.0 * (LMST(MODJD, LAMBDA) - Ra);
  SIN_HSUN := SPHI * SN(Dec) + CPHI * CS(Dec) * CS(TAU);
  DARKNESS := (SIN_HSUN < -0.10);
end;

(* ----------------------------------------------------------------------- *)
(* Examine: *)
(* checks whether an occultation takes place, calculates the *)
(* circumstances and prints the results *)
(* *)
(* T1,T2                : search interval in Julian cent. since J2000 *)
(* RAPOLY,DEPOLY,RPOLY  : Chebyshev approximations of lunar coordinates *)
(* ETDIFUT              : ET-UT in sec *)
(* LAMBDA,RCPHI,RSPHI   : geocentric coordinates of the observer *)
(* RA_STAR,DE_STAR      : star's coordinates *)
(* NAME                 : star's name *)
(* ----------------------------------------------------------------------- *)

procedure Examine(T1, T2: Double; RAPOLY, DEPOLY, RPOLY: TPolynomCheb;
  ETDIFUT, LAMBDA, RCPHI, RSPHI, RA_STAR, DE_STAR: Double; NAME: NAME_STRING);

const
  H_MIN = 5.0; (* minimum altitude above the horizon (deg) *)

var
  DAY, MONTH, YEAR, H, M, I: integer;
  S, HOUR, T_CONJ_ET, MJD_UT_IN, MJD_UT_OUT: Double;
  POS_IN, POS_OUT, H_IN, H_OUT, A_IN, A_OUT, B_IN, B_OUT: Double;
  CONJ, TAKES_PLACE: Boolean;

begin

  (* test for conjunction in Ra and find time of conjunction *)

  Conjunct(T1, T2, RAPOLY, DEPOLY, RPOLY, RA_STAR, DE_STAR, CONJ, T_CONJ_ET);

  if CONJ then

  begin

    (* check a possible occultation for the given observing site *)
    (* and calculate times of contact, altitudes and longitude and *)
    (* latitude coefficients *)

    Shadow(RAPOLY, DEPOLY, RPOLY, T_CONJ_ET, ETDIFUT, LAMBDA, RCPHI, RSPHI, RA_STAR, DE_STAR,
      TAKES_PLACE, MJD_UT_IN, MJD_UT_OUT, POS_IN, POS_OUT, H_IN, H_OUT, A_IN, A_OUT, B_IN, B_OUT);

    (* print results if the occultation takes place during the *)
    (* night and high enough above the horizon *)

    if TAKES_PLACE then
      if ((H_IN > H_MIN) or (H_OUT > H_MIN)) then
        if DARKNESS((MJD_UT_IN + MJD_UT_OUT) / 2.0, LAMBDA, RCPHI, RSPHI) then
        begin

          (* disappearance *)
          CalDat(MJD_UT_IN, DAY, MONTH, YEAR, HOUR);
          DMS(HOUR, H, M, S);
          write((YEAR MOD 100):3, '/', MONTH:2, '/', DAY:2, H:5, M:3, Trunc(S + 0.5):3, '   D  ',
            Trunc(POS_IN + 0.5):5, Trunc(H_IN + 0.5):6, A_IN:8:1, B_IN:6:1, ' ':3);
          for I := 1 to NAME_LENGTH do
            write(NAME[I]);
          writeln;

          (* reapparence *)
          CalDat(MJD_UT_OUT, DAY, MONTH, YEAR, HOUR);
          DMS(HOUR, H, M, S);
          writeln((YEAR MOD 100):3, '/', MONTH:2, '/', DAY:2, H:5, M:3, Trunc(S + 0.5):3, '   R  ',
            Trunc(POS_OUT + 0.5):5, Trunc(H_OUT + 0.5):6, A_OUT:8:1, B_OUT:6:1);

        end;

  end;

end;

(* ----------------------------------------------------------------------- *)

begin (* main program *)

  (* read search interval and geographic coordinates *)

  Get_Input(T_begin, T_END, ETDIFUT, LAMBDA, PHI);

  (* calculate geocentric coordinates of the observer *)

  Site(PHI, RCPHI, RSPHI);

  (* search occultations in subsequent time intervals *)

  T2 := T_begin;

  repeat

    T1 := T2;
    T2 := T1 + T_SEARCH;

    (* approximate lunar coordinates by Chebyshev polynomials *)

    T_Fit_Moon(T1 - TOVLAP, T2 + TOVLAP, MAX_TP_DEG, RAPOLY, DEPOLY, RPOLY);

    (* print header *)

    HEADER;

    (* open star catalogue file, read epoch and equinox *)

    (* Reset ( OCCINP ); *)                           (* Standard Pascal *)
    Assign(OCCINP, 'OCCINP.DAT');
    Reset(OCCINP); (* Turbo Pascal *)

    readln(OCCINP, T_EPOCH, T_EQX);
    T_EQX := (T_EQX - 2000.0) / 100.0;
    T_EPOCH := (T_EPOCH - 2000.0) / 100.0;

    (* calculate transformation matrix between the mean equinox of the *)
    (* star catalog and the true equinox of the search interval center *)

    TM := (T1 + T2) / 2.0;
    PN_Matrix(T_EQX, TM, PNMAT);

    (* heliocentric velocity of the earth for calculation of aberration *)

    Aberrat(TM, VX, VY, VZ);

    (* loop through list of stars and search for possible occultations *)

    while not EOF(OCCINP) do
    begin
      (* read new star coordinates *)
      GETSTAR(T_EPOCH, TM, RA_STAR, DE_STAR, NAME);
      (* calculate apparent coordinates *)
      Apparent(PNMAT, VX, VY, VZ, RA_STAR, DE_STAR);
      (* check for occultation *)
      Examine(T1, T2, RAPOLY, DEPOLY, RPOLY, ETDIFUT, LAMBDA, RCPHI, RSPHI, RA_STAR, DE_STAR, NAME);
    end;
    writeln;

    CLOSE(OCCINP); (* Close input file; Turbo Pascal only *)

  until (T2 >= T_END);

end.
