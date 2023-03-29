unit Apc.Planets;

interface

// Calculation of planetary positions

uses
  System.Math,
  Apc.Kepler,
  Apc.Mathem,
  Apc.Spheric,
  Apc.Sun;


type
  PlanetType = (Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto);

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* PlanetPos:                                                             *)
(*                                                                       *)
(*   Computes the position of a planet assuming Keplerian orbits. Mean   *)
(*   elements at epoch J2000 are used for Mercury to Mars; osculating    *)
(*   elements at epoch 1995/10/10 (JD 2450000.5) are used for Jupiter to *)
(*   Pluto. Relative accuracy approx. 0.001 between 1990 and 2000.       *)
(*                                                                       *)
(*   Planet  Name of the planet                                          *)
(*   T       Time in Julian centuries since J2000                        *)
(*   X,Y,Z   Ecliptic coordinates (equinox J2000)                        *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure PlanetPos(Planet: PlanetType; T: Double; var X, Y, Z: Double);

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* GeoCentric: geocentric coordinates (geometric and light-time corrected)   *)
(*                                                                       *)
(*   T:        time in Julian centuries since J2000                      *)
(*             T=(JD-2451545.0)/36525.0                                  *)
(*   LP,BP,RP: ecliptic heliocentric coordinates of the planet           *)
(*   LS,BS,RS: ecliptic geocentric   coordinates of the sun              *)
(*                                                                       *)
(*   IPLAN:    planet (0=Sun,1=Mercury,2=Venus,3=Earth,...,9=Pluto       *)
(*   IMODE:    desired type of coordinates (see description of X,Y,Z)    *)
(*             (0=geometric,1=astrometric,2=apparent)                    *)
(*   XP,YP,ZP: ecliptic heliocentric coordinates of the planet           *)
(*   XS,YS,ZS: ecliptic geocentric coordinates of the Sun                *)
(*   X, Y, Z : ecliptic geocentric cordinates of the planet (geometric   *)
(*             if IMODE=0, astrometric if IMODE=1, apparent if IMODE=2)  *)
(*   DELTA0:   geocentric distance (geometric)                           *)
(*                                                                       *)
(*   (all angles in degrees, distances in AU)                            *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure GeoCentric(T, LP, BP, RP, LS, BS, RS: Double; IPLAN, IMODE: integer;
  var XP, YP, ZP, XS, YS, ZS, X, Y, Z, DELTA0: Double);

(*-----------------------------------------------------------------------*)
(* MercuryPos: Mercury; ecliptic coordinates L,B,R (in deg and AU)       *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure MercuryPos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* VenusPos: Venus; ecliptic coordinates L,B,R (in deg and AU)           *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure VenusPos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* MarsPos: Mars; ecliptic coordinates L,B,R (in deg and AU)             *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure MarsPos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* JupiterPos: Jupiter; ecliptic coordinates L,B,R (in deg and AU)       *)
(*         equinox of date                                               *)
(*         T: time in Julian centuries since J2000                       *)
(*            = (JED-2451545.0)/36525                                    *)
(*-----------------------------------------------------------------------*)
procedure JupiterPos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* SaturnPos: Saturn; ecliptic coordinates L,B,R (in deg and AU)         *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure SaturnPos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* UranusPos: Uranus; ecliptic coordinates L,B,R (in deg and AU)         *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure UranusPos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* NeptunePos: Neptune; ecliptic coordinates L,B,R (in deg and AU)       *)
(*         equinox of date                                               *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure NeptunePos(T: Double; var L, B, R: Double);

(*-----------------------------------------------------------------------*)
(* PlutoPos: Pluto; ecliptic coordinates L,B,R (in deg and AU)           *)
(*         equinox of date; only valid between 1890 and 2100!!           *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
procedure PlutoPos(T: Double; var L, B, R: Double);

//=====================================================
implementation
//=====================================================


procedure PlanetPos(Planet: PlanetType; T: Double; var X, Y, Z: Double);
var
  A, E, M, O, I, W, N: Double;
  XX, YY, VVX, VVY: Double;
  T0: Double;
  PQR: Double33;
begin
  // Heliocentric ecliptic orbital elements for equinox J2000
  case Planet of
    Mercury:
      begin
        A := 0.387099;
        E := 0.205634;
        M := 174.7947;
        N := 149472.6738;
        O := 48.331;
        I := 7.0048;
        W := 77.4552;
      end;
    Venus:
      begin
        A := 0.723332;
        E := 0.006773;
        M := 50.4071;
        N := 58517.8149;
        O := 76.680;
        I := 3.3946;
        W := 131.5718;
      end;
    Earth:
      begin
        A := 1.000000;
        E := 0.016709;
        M := 357.5256;
        N := 35999.3720;
        O := 174.876;
        I := 0.0000;
        W := 102.9400;
      end;
    Mars:
      begin
        A := 1.523692;
        E := 0.093405;
        M := 19.3879;
        N := 19140.3023;
        O := 49.557;
        I := 1.8496;
        W := 336.0590;
      end;
    Jupiter:
      begin
        A := 5.202437;
        E := 0.048402;
        M := 250.3274;
        N := 3035.2275;
        O := 100.4683;
        I := 1.3047;
        W := 15.7192;
      end;
    Saturn:
      begin
        A := 9.551712;
        E := 0.052340;
        M := 267.2465;
        N := 1219.6465;
        O := 113.6439;
        I := 2.4855;
        W := 90.9682;
      end;
    Uranus:
      begin
        A := 19.293108;
        E := 0.044846;
        M := 118.4320;
        N := 424.8150;
        O := 74.0903;
        I := 0.7733;
        W := 176.6152;
      end;
    Neptune:
      begin
        A := 30.257162;
        E := 0.007985;
        M := 292.4716;
        N := 216.3047;
        O := 131.7750;
        I := 1.7700;
        W := 3.0962;
      end;
    Pluto:
      begin
        A := 39.783607;
        E := 0.254351;
        M := 8.2304;
        N := 143.4629;
        O := 110.3865;
        I := 17.1201;
        W := 224.7424;
      end;
  end;
  case Planet of
    Mercury, Venus, Earth, Mars:
      T0 := 0.0;
    Jupiter, Saturn, Uranus, Neptune, Pluto:
      T0 := -0.042286105;
  end;
  M := M + N * (T - T0);

  // Cartesian coordinates mean ecliptic and equinox J2000
  GaussVec(O, I, W - O, PQR);
  Ellip(M, A, E, XX, YY, VVX, VVY);
  Orb2Ecl(XX, YY, PQR, X, Y, Z);
end; // Position

(* --------------------------------------------------------------------- *)

procedure GeoCentric(T, LP, BP, RP, LS, BS, RS: Double; IPLAN, IMODE: integer;
  var XP, YP, ZP, XS, YS, ZS, X, Y, Z, DELTA0: Double);

const
  P2 = 6.283185307;
var
  DL, DB, DR, DLS, DBS, DRS, FAC: Double;
  VX, VY, VZ, VXS, VYS, VZS, M: Double;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1;
    Frac := X
  end;

  (*sub*)procedure PosVel(L, B, R, DL, DB, DR: Double; var X, Y, Z, VX, VY, VZ: Double);
  var
    CL, SL, CB, SB: Double;
  begin
    CL := CS(L);
    SL := SN(L);
    CB := CS(B);
    SB := SN(B);
    X := R * CL * CB;
    VX := DR * CL * CB - DL * R * SL * CB - DB * R * CL * SB;
    Y := R * SL * CB;
    VY := DR * SL * CB + DL * R * CL * CB - DB * R * SL * SB;
    Z := R * SB;
    VZ := DR * SB + DB * R * CB;
  end;

begin
  DL := 0.0;
  DB := 0.0;
  DR := 0.0;
  DLS := 0.0;
  DBS := 0.0;
  DRS := 0.0;

  if (IMODE > 0) then
  begin
    M := P2 * Frac(0.9931266 + 99.9973604 * T); // Sun
    DLS := 172.00 + 5.75 * Sin(M);
    DRS := 2.87 * Cos(M);
    DBS := 0.0;

    // dl,db in 1e-4 rad/d, dr in 1e-4 AU/d
    case IPLAN of
      0: // Sun
        begin
          DL := 0.0;
          DB := 0.0;
          DR := 0.0;
        end;
      1: // Mercury
        begin
          M := P2 * Frac(0.4855407 + 415.2014314 * T);
          DL := 714.00 + 292.66 * Cos(M) + 71.96 * Cos(2 * M) + 18.16 * Cos(3 * M) + 4.61 *
            Cos(4 * M) + 3.81 * Sin(2 * M) + 2.43 * Sin(3 * M) + 1.08 * Sin(4 * M);
          DR := 55.94 * Sin(M) + 11.36 * Sin(2 * M) + 2.60 * Sin(3 * M);
          DB := 73.40 * Cos(M) + 29.82 * Cos(2 * M) + 10.22 * Cos(3 * M) + 3.28 * Cos(4 * M) - 40.44
            * Sin(M) - 16.55 * Sin(2 * M) - 5.56 * Sin(3 * M) - 1.72 * Sin(4 * M);
        end;
      2: // Venus
        begin
          M := P2 * Frac(0.1400197 + 162.5494552 * T);
          DL := 280.00 + 3.79 * Cos(M);
          DR := 1.37 * Sin(M);
          DB := 9.54 * Cos(M) - 13.57 * Sin(M);
        end;
      3: // Earth
        begin
          DL := DLS;
          DR := DRS;
          DB := -DBS;
        end;
      4: // Mars
        begin
          M := P2 * Frac(0.0538553 + 53.1662736 * T);
          DL := 91.50 + 17.07 * Cos(M) + 2.03 * Cos(2 * M);
          DR := 12.98 * Sin(M) + 1.21 * Cos(2 * M);
          DB := 0.83 * Cos(M) + 2.80 * Sin(M);
        end;
      5: // Jupiter
        begin
          M := P2 * Frac(0.0565314 + 8.4302963 * T);
          DL := 14.50 + 1.41 * Cos(M);
          DR := 3.66 * Sin(M);
          DB := 0.33 * Sin(M);
        end;
      6: // Saturn
        begin
          M := P2 * Frac(0.8829867 + 3.3947688 * T);
          DL := 5.84 + 0.65 * Cos(M);
          DR := 3.09 * Sin(M);
          DB := 0.24 * Cos(M);
        end;
      7: // Uranus
        begin
          M := P2 * Frac(0.3967117 + 1.1902849 * T);
          DL := 2.05 + 0.19 * Cos(M);
          DR := 1.86 * Sin(M);
          DB := -0.03 * Sin(M);
        end;
      8: // Neptune
        begin
          M := P2 * Frac(0.7214906 + 0.6068526 * T);
          DL := 1.04 + 0.02 * Cos(M);
          DR := 0.27 * Sin(M);
          DB := 0.03 * Sin(M);
        end;
      9: // Pluto
        begin
          M := P2 * Frac(0.0385795 + 0.4026667 * T);
          DL := 0.69 + 0.34 * Cos(M) + 0.12 * Cos(2 * M) + 0.05 * Cos(3 * M);
          DR := 6.66 * Sin(M) + 1.64 * Sin(2 * M);
          DB := -0.08 * Cos(M) - 0.17 * Sin(M) - 0.09 * Sin(2 * M);
        end;
    end;
  end;

  PosVel(LS, BS, RS, DLS, DBS, DRS, XS, YS, ZS, VXS, VYS, VZS);
  PosVel(LP, BP, RP, DL, DB, DR, XP, YP, ZP, VX, VY, VZ);
  X := XP + XS;
  Y := YP + YS;
  Z := ZP + ZS;
  DELTA0 := Sqrt(X * X + Y * Y + Z * Z);
  if IPLAN = 3 then
  begin
    X := 0.0;
    Y := 0.0;
    Z := 0.0;
    DELTA0 := 0.0
  end;

  FAC := 0.00578 * DELTA0 * 1E-4;
  case IMODE of
    1:
      begin
        X := X - FAC * VX;
        Y := Y - FAC * VY;
        Z := Z - FAC * VZ;
      end;
    2:
      begin
        X := X - FAC * (VX + VXS);
        Y := Y - FAC * (VY + VYS);
        Z := Z - FAC * (VZ + VZS);
      end;
  end;
end;


(* ----------------------------------------------------------------------- *)

procedure MercuryPos(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C1, S1: array [-1 .. 9] of Double;
  C, S: array [-5 .. 0] of Double;
  M1, M2, M3, M5, M6: Double;
  U, V, DL, DR, DB: Double;
  I: integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I1, I, IT: integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C1[I1], S1[I1], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertVenus; // Kepler terms and perturbations by Venus
  var
    I: integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M2);
    S[-1] := -Sin(M2);
    for I := -1 downto -4 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, 0, 0, 259.74, 84547.39, -78342.34, 0.01, 11683.22, 21203.79);
    Term(1, 0, 1, 2.30, 5.04, -7.52, 0.02, 138.55, -71.01);
    Term(1, 0, 2, 0.01, -0.01, 0.01, 0.01, -0.19, -0.54);
    Term(2, 0, 0, -549.71, 10394.44, -7955.45, 0.00, 2390.29, 4306.79);
    Term(2, 0, 1, -4.77, 8.97, -1.53, 0.00, 28.49, -14.18);
    Term(2, 0, 2, 0.00, 0.00, 0.00, 0.00, -0.04, -0.11);
    Term(3, 0, 0, -234.04, 1748.74, -1212.86, 0.00, 535.41, 984.33);
    Term(3, 0, 1, -2.03, 3.48, -0.35, 0.00, 6.56, -2.91);
    Term(4, 0, 0, -77.64, 332.63, -219.23, 0.00, 124.40, 237.03);
    Term(4, 0, 1, -0.70, 1.10, -0.08, 0.00, 1.59, -0.59);
    Term(5, 0, 0, -23.59, 67.28, -43.54, 0.00, 29.44, 58.77);
    Term(5, 0, 1, -0.23, 0.32, -0.02, 0.00, 0.39, -0.11);
    Term(6, 0, 0, -6.86, 14.06, -9.18, 0.00, 7.03, 14.84);
    Term(6, 0, 1, -0.07, 0.09, -0.01, 0.00, 0.10, -0.02);
    Term(7, 0, 0, -1.94, 2.98, -2.02, 0.00, 1.69, 3.80);
    Term(8, 0, 0, -0.54, 0.63, -0.46, 0.00, 0.41, 0.98);
    Term(9, 0, 0, -0.15, 0.13, -0.11, 0.00, 0.10, 0.25);
    Term(-1, -2, 0, -0.17, -0.06, -0.05, 0.14, -0.06, -0.07);
    Term(0, -1, 0, 0.24, -0.16, -0.11, -0.16, 0.04, -0.01);
    Term(0, -2, 0, -0.68, -0.25, -0.26, 0.73, -0.16, -0.18);
    Term(0, -5, 0, 0.37, 0.08, 0.06, -0.28, 0.13, 0.12);
    Term(1, -1, 0, 0.58, -0.41, 0.26, 0.36, 0.01, -0.01);
    Term(1, -2, 0, -3.51, -1.23, 0.23, -0.63, -0.05, -0.06);
    Term(1, -3, 0, 0.08, 0.53, -0.11, 0.04, 0.02, -0.09);
    Term(1, -5, 0, 1.44, 0.31, 0.30, -1.39, 0.34, 0.29);
    Term(2, -1, 0, 0.15, -0.11, 0.09, 0.12, 0.02, -0.04);
    Term(2, -2, 0, -1.99, -0.68, 0.65, -1.91, -0.20, 0.03);
    Term(2, -3, 0, -0.34, -1.28, 0.97, -0.26, 0.03, 0.03);
    Term(2, -4, 0, -0.33, 0.35, -0.13, -0.13, -0.01, 0.00);
    Term(2, -5, 0, 7.19, 1.56, -0.05, 0.12, 0.06, 0.05);
    Term(3, -2, 0, -0.52, -0.18, 0.13, -0.39, -0.16, 0.03);
    Term(3, -3, 0, -0.11, -0.42, 0.36, -0.10, -0.05, -0.05);
    Term(3, -4, 0, -0.19, 0.22, -0.23, -0.20, -0.01, 0.02);
    Term(3, -5, 0, 2.77, 0.49, -0.45, 2.56, 0.40, -0.12);
    Term(4, -5, 0, 0.67, 0.12, -0.09, 0.47, 0.24, -0.08);
    Term(5, -5, 0, 0.18, 0.03, -0.02, 0.12, 0.09, -0.03);
  end;

  (*sub*)procedure PertEarth; // perturbations by the Earth
  var
    I: integer;
  begin
    C[-1] := Cos(M3);
    S[-1] := -Sin(M3);
    for I := -1 downto -3 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(0, -4, 0, -0.11, -0.07, -0.08, 0.11, -0.02, -0.04);
    Term(1, -1, 0, 0.10, -0.20, 0.15, 0.07, 0.00, 0.00);
    Term(1, -2, 0, -0.35, 0.28, -0.13, -0.17, -0.01, 0.00);
    Term(1, -4, 0, -0.67, -0.45, 0.00, 0.01, -0.01, -0.01);
    Term(2, -2, 0, -0.20, 0.16, -0.16, -0.20, -0.01, 0.02);
    Term(2, -3, 0, 0.13, -0.02, 0.02, 0.14, 0.01, 0.00);
    Term(2, -4, 0, -0.33, -0.18, 0.17, -0.31, -0.04, 0.00);
  end;

  (*sub*)procedure PertJupiter; // perturbations by Jupiter
  var
    I: integer;
  begin
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    for I := -1 downto -2 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(-1, -1, 0, -0.08, 0.16, 0.15, 0.08, -0.04, 0.01);
    Term(-1, -2, 0, 0.10, -0.06, -0.07, -0.12, 0.07, -0.01);
    Term(0, -1, 0, -0.31, 0.48, -0.02, 0.13, -0.03, -0.02);
    Term(0, -2, 0, 0.42, -0.26, -0.38, -0.50, 0.20, -0.03);
    Term(1, -1, 0, -0.70, 0.01, -0.02, -0.63, 0.00, 0.03);
    Term(1, -2, 0, 2.61, -1.97, 1.74, 2.32, 0.01, 0.01);
    Term(1, -3, 0, 0.32, -0.15, 0.13, 0.28, 0.00, 0.00);
    Term(2, -1, 0, -0.18, 0.01, 0.00, -0.13, -0.03, 0.03);
    Term(2, -2, 0, 0.75, -0.56, 0.45, 0.60, 0.08, -0.17);
    Term(3, -2, 0, 0.20, -0.15, 0.10, 0.14, 0.04, -0.08);
  end;

  (*sub*)procedure PertSaturn; // perturbations by Saturn
  begin
    C[-2] := Cos(2 * M6);
    S[-2] := -Sin(2 * M6);
    Term(1, -2, 0, -0.19, 0.33, 0.00, 0.00, 0.00, 0.00);
  end;

begin // MercuryPos

  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M1 := P2 * Frac(0.4855407 + 415.2014314 * T);
  M2 := P2 * Frac(0.1394222 + 162.5490444 * T);
  M3 := P2 * Frac(0.9937861 + 99.9978139 * T);
  M5 := P2 * Frac(0.0558417 + 8.4298417 * T);
  M6 := P2 * Frac(0.8823333 + 3.3943333 * T);
  C1[0] := 1.0;
  S1[0] := 0.0;
  C1[1] := Cos(M1);
  S1[1] := Sin(M1);
  C1[-1] := C1[1];
  S1[-1] := -S1[1];
  for I := 2 to 9 do
    AddThe(C1[I - 1], S1[I - 1], C1[1], S1[1], C1[I], S1[I]);
  PertVenus;
  PertEarth;
  PertJupiter;
  PertSaturn;
  DL := DL + (2.8 + 3.2 * T);
  L := 360.0 * Frac(0.2151379 + M1 / P2 + ((5601.7 + 1.1 * T) * T + DL) / 1296.0E3);
  R := 0.3952829 + 0.0000016 * T + DR * 1.0E-6;
  B := (-2522.15 + (-30.18 + 0.04 * T) * T + DB) / 3600.0;
end; // MercuryPos

(* ----------------------------------------------------------------------- *)

procedure VenusPos(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C2, S2: array [0 .. 8] of Double;
  C, S: array [-8 .. 0] of Double;
  M1, M2, M3, M4, M5, M6: Double;
  U, V, DL, DR, DB: Double;
  I: integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I1, I, IT: integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C2[I1], S2[I1], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertMercury; // perturbations by Mercury
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M1);
    S[-1] := -Sin(M1);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(1, -1, 0, 0.00, 0.00, 0.06, -0.09, 0.01, 0.00);
    Term(2, -1, 0, 0.25, -0.09, -0.09, -0.27, 0.00, 0.00);
    Term(4, -2, 0, -0.07, -0.08, -0.14, 0.14, -0.01, -0.01);
    Term(5, -2, 0, -0.35, 0.08, 0.02, 0.09, 0.00, 0.00);
  end;

  (*sub*)procedure PertEarth; // Kepler terms and perturbations by the Earth
  var
    I: integer;
  begin
    C[-1] := Cos(M3);
    S[-1] := -Sin(M3);
    for I := -1 downto -7 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, 0, 0, 2.37, 2793.23, -4899.07, 0.11, 9995.27, 7027.22);
    Term(1, 0, 1, 0.10, -19.65, 34.40, 0.22, 64.95, -86.10);
    Term(1, 0, 2, 0.06, 0.04, -0.07, 0.11, -0.55, -0.07);
    Term(2, 0, 0, -170.42, 73.13, -16.59, 0.00, 67.71, 47.56);
    Term(2, 0, 1, 0.93, 2.91, 0.23, 0.00, -0.03, -0.92);
    Term(3, 0, 0, -2.31, 0.90, -0.08, 0.00, 0.04, 2.09);
    Term(1, -1, 0, -2.38, -4.27, 3.27, -1.82, 0.00, 0.00);
    Term(1, -2, 0, 0.09, 0.00, -0.08, 0.05, -0.02, -0.25);
    Term(2, -2, 0, -9.57, -5.93, 8.57, -13.83, -0.01, -0.01);
    Term(2, -3, 0, -2.47, -2.40, 0.83, -0.95, 0.16, 0.24);
    Term(3, -2, 0, -0.09, -0.05, 0.08, -0.13, -0.28, 0.12);
    Term(3, -3, 0, 7.12, 0.32, -0.62, 13.76, -0.07, 0.01);
    Term(3, -4, 0, -0.65, -0.17, 0.18, -0.73, 0.10, 0.05);
    Term(3, -5, 0, -1.08, -0.95, -0.17, 0.22, -0.03, -0.03);
    Term(4, -3, 0, 0.06, 0.00, -0.01, 0.08, 0.14, -0.18);
    Term(4, -4, 0, 0.93, -0.46, 1.06, 2.13, -0.01, 0.01);
    Term(4, -5, 0, -1.53, 0.38, -0.64, -2.54, 0.27, 0.00);
    Term(4, -6, 0, -0.17, -0.05, 0.03, -0.11, 0.02, 0.00);
    Term(5, -5, 0, 0.18, -0.28, 0.71, 0.47, -0.02, 0.04);
    Term(5, -6, 0, 0.15, -0.14, 0.30, 0.31, -0.04, 0.03);
    Term(5, -7, 0, -0.08, 0.02, -0.03, -0.11, 0.01, 0.00);
    Term(5, -8, 0, -0.23, 0.00, 0.01, -0.04, 0.00, 0.00);
    Term(6, -6, 0, 0.01, -0.14, 0.39, 0.04, 0.00, -0.01);
    Term(6, -7, 0, 0.02, -0.05, 0.12, 0.04, -0.01, 0.01);
    Term(6, -8, 0, 0.10, -0.10, 0.19, 0.19, -0.02, 0.02);
    Term(7, -7, 0, -0.03, -0.06, 0.18, -0.08, 0.00, 0.00);
    Term(8, -8, 0, -0.03, -0.02, 0.06, -0.08, 0.00, 0.00);
  end;

  (*sub*)procedure PertMars; // perturbations by Mars
  var
    I: integer;
  begin
    C[-1] := Cos(M4);
    S[-1] := -Sin(M4);
    for I := -1 downto -2 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, -3, 0, -0.65, 1.02, -0.04, -0.02, -0.02, 0.00);
    Term(2, -2, 0, -0.05, 0.04, -0.09, -0.10, 0.00, 0.00);
    Term(2, -3, 0, -0.50, 0.45, -0.79, -0.89, 0.01, 0.03);
  end;

  (*sub*)procedure PertJupiter; // perturbations by Jupiter
  var
    I: integer;
  begin
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    for I := -1 downto -2 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(0, -1, 0, -0.05, 1.56, 0.16, 0.04, -0.08, -0.04);
    Term(1, -1, 0, -2.62, 1.40, -2.35, -4.40, 0.02, 0.03);
    Term(1, -2, 0, -0.47, -0.08, 0.12, -0.76, 0.04, -0.18);
    Term(2, -2, 0, -0.73, -0.51, 1.27, -1.82, -0.01, 0.01);
    Term(2, -3, 0, -0.14, -0.10, 0.25, -0.34, 0.00, 0.00);
    Term(3, -3, 0, -0.01, 0.04, -0.11, -0.02, 0.00, 0.00);
  end;

  (*sub*)procedure PertSaturn; // perturbations by Saturn
  begin
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    Term(0, -1, 0, 0.00, 0.21, 0.00, 0.00, 0.00, -0.01);
    Term(1, -1, 0, -0.11, -0.14, 0.24, -0.20, 0.01, 0.00);
  end;

begin // VenusPos
  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M1 := P2 * Frac(0.4861431 + 415.2018375 * T);
  M2 := P2 * Frac(0.1400197 + 162.5494552 * T);
  M3 := P2 * Frac(0.9944153 + 99.9982208 * T);
  M4 := P2 * Frac(0.0556297 + 53.1674631 * T);
  M5 := P2 * Frac(0.0567028 + 8.4305083 * T);
  M6 := P2 * Frac(0.8830539 + 3.3947206 * T);
  C2[0] := 1.0;
  S2[0] := 0.0;
  C2[1] := Cos(M2);
  S2[1] := Sin(M2);
  for I := 2 to 8 do
    AddThe(C2[I - 1], S2[I - 1], C2[1], S2[1], C2[I], S2[I]);
  PertMercury;
  PertEarth;
  PertMars;
  PertJupiter;
  PertSaturn;
  DL := DL + 2.74 * Sin(P2 * (0.0764 + 0.4174 * T)) + 0.27 * Sin(P2 * (0.9201 + 0.3307 * T));
  DL := DL + (1.9 + 1.8 * T);
  L := 360.0 * Frac(0.3654783 + M2 / P2 + ((5071.2 + 1.1 * T) * T + DL) / 1296.0E3);
  R := 0.7233482 - 0.0000002 * T + DR * 1.0E-6;
  B := (-67.70 + (0.04 + 0.01 * T) * T + DB) / 3600.0;

end; // VenusPos


//-----------------------------------------------------------------------

procedure MarsPos(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C4, S4: array [-2 .. 16] of Double;
  C, S: array [-9 .. 0] of Double;
  M2, M3, M4, M5, M6: Double;
  U, V, DL, DR, DB: Double;
  I: integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I1, I, IT: integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C4[I1], S4[I1], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertVenus; // perturbations by Venus
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M2);
    S[-1] := -Sin(M2);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(0, -1, 0, -0.01, -0.03, 0.10, -0.04, 0.00, 0.00);
    Term(1, -1, 0, 0.05, 0.10, -2.08, 0.75, 0.00, 0.00);
    Term(2, -1, 0, -0.25, -0.57, -2.58, 1.18, 0.05, -0.04);
    Term(2, -2, 0, 0.02, 0.02, 0.13, -0.14, 0.00, 0.00);
    Term(3, -1, 0, 3.41, 5.38, 1.87, -1.15, 0.01, -0.01);
    Term(3, -2, 0, 0.02, 0.02, 0.11, -0.13, 0.00, 0.00);
    Term(4, -1, 0, 0.32, 0.49, -1.88, 1.21, -0.07, 0.07);
    Term(4, -2, 0, 0.03, 0.03, 0.12, -0.14, 0.00, 0.00);
    Term(5, -1, 0, 0.04, 0.06, -0.17, 0.11, -0.01, 0.01);
    Term(5, -2, 0, 0.11, 0.09, 0.35, -0.43, -0.01, 0.01);
    Term(6, -2, 0, -0.36, -0.28, -0.20, 0.25, 0.00, 0.00);
    Term(7, -2, 0, -0.03, -0.03, 0.11, -0.13, 0.00, -0.01);
  end;

  (*sub*)procedure PertEarth; // Kepler terms and perturbations by the Earth
  var
    I: integer;
  begin
    C[-1] := Cos(M3);
    S[-1] := -Sin(M3);
    for I := -1 downto -8 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(1, 0, 0, -5.32, 38481.97, -141856.04, 0.40, -6321.67, 1876.89);
    Term(1, 0, 1, -1.12, 37.98, -138.67, -2.93, 37.28, 117.48);
    Term(1, 0, 2, -0.32, -0.03, 0.12, -1.19, 1.04, -0.40);
    Term(2, 0, 0, 28.28, 2285.80, -6608.37, 0.00, -589.35, 174.81);
    Term(2, 0, 1, 1.64, 3.37, -12.93, 0.00, 2.89, 11.10);
    Term(2, 0, 2, 0.00, 0.00, 0.00, 0.00, 0.10, -0.03);
    Term(3, 0, 0, 5.31, 189.29, -461.81, 0.00, -61.98, 18.53);
    Term(3, 0, 1, 0.31, 0.35, -1.36, 0.00, 0.25, 1.19);
    Term(4, 0, 0, 0.81, 17.96, -38.26, 0.00, -6.88, 2.08);
    Term(4, 0, 1, 0.05, 0.04, -0.15, 0.00, 0.02, 0.14);
    Term(5, 0, 0, 0.11, 1.83, -3.48, 0.00, -0.79, 0.24);
    Term(6, 0, 0, 0.02, 0.20, -0.34, 0.00, -0.09, 0.03);
    Term(-1, -1, 0, 0.09, 0.06, 0.14, -0.22, 0.02, -0.02);
    Term(0, -1, 0, 0.72, 0.49, 1.55, -2.31, 0.12, -0.10);
    Term(1, -1, 0, 7.00, 4.92, 13.93, -20.48, 0.08, -0.13);
    Term(2, -1, 0, 13.08, 4.89, -4.53, 10.01, -0.05, 0.13);
    Term(2, -2, 0, 0.14, 0.05, -0.48, -2.66, 0.01, 0.14);
    Term(3, -1, 0, 1.38, 0.56, -2.00, 4.85, -0.01, 0.19);
    Term(3, -2, 0, -6.85, 2.68, 8.38, 21.42, 0.00, 0.03);
    Term(3, -3, 0, -0.08, 0.20, 1.20, 0.46, 0.00, 0.00);
    Term(4, -1, 0, 0.16, 0.07, -0.19, 0.47, -0.01, 0.05);
    Term(4, -2, 0, -4.41, 2.14, -3.33, -7.21, -0.07, -0.09);
    Term(4, -3, 0, -0.12, 0.33, 2.22, 0.72, -0.03, -0.02);
    Term(4, -4, 0, -0.04, -0.06, -0.36, 0.23, 0.00, 0.00);
    Term(5, -2, 0, -0.44, 0.21, -0.70, -1.46, -0.06, -0.07);
    Term(5, -3, 0, 0.48, -2.60, -7.25, -1.37, 0.00, 0.00);
    Term(5, -4, 0, -0.09, -0.12, -0.66, 0.50, 0.00, 0.00);
    Term(5, -5, 0, 0.03, 0.00, 0.01, -0.17, 0.00, 0.00);
    Term(6, -2, 0, -0.05, 0.03, -0.07, -0.15, -0.01, -0.01);
    Term(6, -3, 0, 0.10, -0.96, 2.36, 0.30, 0.04, 0.00);
    Term(6, -4, 0, -0.17, -0.20, -1.09, 0.94, 0.02, -0.02);
    Term(6, -5, 0, 0.05, 0.00, 0.00, -0.30, 0.00, 0.00);
    Term(7, -3, 0, 0.01, -0.10, 0.32, 0.04, 0.02, 0.00);
    Term(7, -4, 0, 0.86, 0.77, 1.86, -2.01, 0.01, -0.01);
    Term(7, -5, 0, 0.09, -0.01, -0.05, -0.44, 0.00, 0.00);
    Term(7, -6, 0, -0.01, 0.02, 0.10, 0.08, 0.00, 0.00);
    Term(8, -4, 0, 0.20, 0.16, -0.53, 0.64, -0.01, 0.02);
    Term(8, -5, 0, 0.17, -0.03, -0.14, -0.84, 0.00, 0.01);
    Term(8, -6, 0, -0.02, 0.03, 0.16, 0.09, 0.00, 0.00);
    Term(9, -5, 0, -0.55, 0.15, 0.30, 1.10, 0.00, 0.00);
    Term(9, -6, 0, -0.02, 0.04, 0.20, 0.10, 0.00, 0.00);
    Term(10, -5, 0, -0.09, 0.03, -0.10, -0.33, 0.00, -0.01);
    Term(10, -6, 0, -0.05, 0.11, 0.48, 0.21, -0.01, 0.00);
    Term(11, -6, 0, 0.10, -0.35, -0.52, -0.15, 0.00, 0.00);
    Term(11, -7, 0, -0.01, -0.02, -0.10, 0.07, 0.00, 0.00);
    Term(12, -6, 0, 0.01, -0.04, 0.18, 0.04, 0.01, 0.00);
    Term(12, -7, 0, -0.05, -0.07, -0.29, 0.20, 0.01, 0.00);
    Term(13, -7, 0, 0.23, 0.27, 0.25, -0.21, 0.00, 0.00);
    Term(14, -7, 0, 0.02, 0.03, -0.10, 0.09, 0.00, 0.00);
    Term(14, -8, 0, 0.05, 0.01, 0.03, -0.23, 0.00, 0.03);
    Term(15, -8, 0, -1.53, 0.27, 0.06, 0.42, 0.00, 0.00);
    Term(16, -8, 0, -0.14, 0.02, -0.10, -0.55, -0.01, -0.02);
    Term(16, -9, 0, 0.03, -0.06, -0.25, -0.11, 0.00, 0.00);
  end;

  (*sub*)procedure PertJupiter; // perturbations by Jupiter
  var
    I: integer;
  begin
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    for I := -1 downto -4 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(-2, -1, 0, 0.05, 0.03, 0.08, -0.14, 0.01, -0.01);
    Term(-1, -1, 0, 0.39, 0.27, 0.92, -1.50, -0.03, -0.06);
    Term(-1, -2, 0, -0.16, 0.03, 0.13, 0.67, -0.01, 0.06);
    Term(-1, -3, 0, -0.02, 0.01, 0.05, 0.09, 0.00, 0.01);
    Term(0, -1, 0, 3.56, 1.13, -5.41, -7.18, -0.25, -0.24);
    Term(0, -2, 0, -1.44, 0.25, 1.24, 7.96, 0.02, 0.31);
    Term(0, -3, 0, -0.21, 0.11, 0.55, 1.04, 0.01, 0.05);
    Term(0, -4, 0, -0.02, 0.02, 0.11, 0.11, 0.00, 0.01);
    Term(1, -1, 0, 16.67, -19.15, 61.00, 53.36, -0.06, -0.07);
    Term(1, -2, 0, -21.64, 3.18, -7.77, -54.64, -0.31, 0.50);
    Term(1, -3, 0, -2.82, 1.45, -2.53, -5.73, 0.01, 0.07);
    Term(1, -4, 0, -0.31, 0.28, -0.34, -0.51, 0.00, 0.00);
    Term(2, -1, 0, 2.15, -2.29, 7.04, 6.94, 0.33, 0.19);
    Term(2, -2, 0, -15.69, 3.31, -15.70, -73.17, -0.17, -0.25);
    Term(2, -3, 0, -1.73, 1.95, -9.19, -7.20, 0.02, -0.03);
    Term(2, -4, 0, -0.01, 0.33, -1.42, 0.08, 0.01, -0.01);
    Term(2, -5, 0, 0.03, 0.03, -0.13, 0.12, 0.00, 0.00);
    Term(3, -1, 0, 0.26, -0.28, 0.73, 0.71, 0.08, 0.04);
    Term(3, -2, 0, -2.06, 0.46, -1.61, -6.72, -0.13, -0.25);
    Term(3, -3, 0, -1.28, -0.27, 2.21, -6.90, -0.04, -0.02);
    Term(3, -4, 0, -0.22, 0.08, -0.44, -1.25, 0.00, 0.01);
    Term(3, -5, 0, -0.02, 0.03, -0.15, -0.08, 0.00, 0.00);
    Term(4, -1, 0, 0.03, -0.03, 0.08, 0.08, 0.01, 0.01);
    Term(4, -2, 0, -0.26, 0.06, -0.17, -0.70, -0.03, -0.05);
    Term(4, -3, 0, -0.20, -0.05, 0.22, -0.79, -0.01, -0.02);
    Term(4, -4, 0, -0.11, -0.14, 0.93, -0.60, 0.00, 0.00);
    Term(4, -5, 0, -0.04, -0.02, 0.09, -0.23, 0.00, 0.00);
    Term(5, -4, 0, -0.02, -0.03, 0.13, -0.09, 0.00, 0.00);
    Term(5, -5, 0, 0.00, -0.03, 0.21, 0.01, 0.00, 0.00);
  end;

  (*sub*)procedure PertSaturn; // perturbations by Saturn
  var
    I: integer;
  begin
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    for I := -1 downto -3 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(-1, -1, 0, 0.03, 0.13, 0.48, -0.13, 0.02, 0.00);
    Term(0, -1, 0, 0.27, 0.84, 0.40, -0.43, 0.01, -0.01);
    Term(0, -2, 0, 0.12, -0.04, -0.33, -0.55, -0.01, -0.02);
    Term(0, -3, 0, 0.02, -0.01, -0.07, -0.08, 0.00, 0.00);
    Term(1, -1, 0, 1.12, 0.76, -2.66, 3.91, -0.01, 0.01);
    Term(1, -2, 0, 1.49, -0.95, 3.07, 4.83, 0.04, -0.05);
    Term(1, -3, 0, 0.21, -0.18, 0.55, 0.64, 0.00, 0.00);
    Term(2, -1, 0, 0.12, 0.10, -0.29, 0.34, -0.01, 0.02);
    Term(2, -2, 0, 0.51, -0.36, 1.61, 2.25, 0.03, 0.01);
    Term(2, -3, 0, 0.10, -0.10, 0.50, 0.43, 0.00, 0.00);
    Term(2, -4, 0, 0.01, -0.02, 0.11, 0.05, 0.00, 0.00);
    Term(3, -2, 0, 0.07, -0.05, 0.16, 0.22, 0.01, 0.01);
  end;

begin // MarsPos

  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M2 := P2 * Frac(0.1382208 + 162.5482542 * T);
  M3 := P2 * Frac(0.9926208 + 99.9970236 * T);
  M4 := P2 * Frac(0.0538553 + 53.1662736 * T);
  M5 := P2 * Frac(0.0548944 + 8.4290611 * T);
  M6 := P2 * Frac(0.8811167 + 3.3935250 * T);
  C4[0] := 1.0;
  S4[0] := 0.0;
  C4[1] := Cos(M4);
  S4[1] := Sin(M4);
  for I := 2 to 16 do
    AddThe(C4[I - 1], S4[I - 1], C4[1], S4[1], C4[I], S4[I]);
  for I := -2 to -1 do
  begin
    C4[I] := C4[-I];
    S4[I] := -S4[-I]
  end;
  PertVenus;
  PertEarth;
  PertJupiter;
  PertSaturn;
  DL := DL + 52.49 * Sin(P2 * (0.1868 + 0.0549 * T)) + 0.61 * Sin(P2 * (0.9220 + 0.3307 * T)) + 0.32
    * Sin(P2 * (0.4731 + 2.1485 * T)) + 0.28 * Sin(P2 * (0.9467 + 0.1133 * T));
  DL := DL + (0.14 + 0.87 * T - 0.11 * T * T);
  L := 360.0 * Frac(0.9334591 + M4 / P2 + ((6615.5 + 1.1 * T) * T + DL) / 1296.0E3);
  R := 1.5303352 + 0.0000131 * T + DR * 1.0E-6;
  B := (596.32 + (-2.92 - 0.10 * T) * T + DB) / 3600.0;
end;

//-----------------------------------------------------------------------

procedure JupiterPos(T: Double; var L, B, R: Double);
const
  P2 = 6.283185307;
var
  C5, S5: array [-1 .. 5] of Double;
  C, S: array [-10 .. 0] of Double;
  M5, M6, M7: Double;
  U, V, DL, DR, DB: Double;
  I: integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I5, I, IT: integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C5[I5], S5[I5], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertSaturn; // Kepler terms and perturbations by Saturn
  var
    I: integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    for I := -1 downto -9 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
    Term(-1, -1, 0, -0.2, 1.4, 2.0, 0.6, 0.1, -0.2);
    Term(0, -1, 0, 9.4, 8.9, 3.9, -8.3, -0.4, -1.4);
    Term(0, -2, 0, 5.6, -3.0, -5.4, -5.7, -2.0, 0.0);
    Term(0, -3, 0, -4.0, -0.1, 0.0, 5.5, 0.0, 0.0);
    Term(0, -5, 0, 3.3, -1.6, -1.6, -3.1, -0.5, -1.2);
    Term(1, 0, 0, -113.1, 19998.6, -25208.2, -142.2, -4670.7, 288.9);
    Term(1, 0, 1, -76.1, 66.9, -84.2, -95.8, 21.6, 29.4);
    Term(1, 0, 2, -0.5, -0.3, 0.4, -0.7, 0.1, -0.1);
    Term(1, -1, 0, 78.8, -14.5, 11.5, 64.4, -0.2, 0.2);
    Term(1, -2, 0, -2.0, -132.4, 28.8, 4.3, -1.7, 0.4);
    Term(1, -2, 1, -1.1, -0.7, 0.2, -0.3, 0.0, 0.0);
    Term(1, -3, 0, -7.5, -6.8, -0.4, -1.1, 0.6, -0.9);
    Term(1, -4, 0, 0.7, 0.7, 0.6, -1.1, 0.0, -0.2);
    Term(1, -5, 0, 51.5, -26.0, -32.5, -64.4, -4.9, -12.4);
    Term(1, -5, 1, -1.2, -2.2, -2.7, 1.5, -0.4, 0.3);
    Term(2, 0, 0, -3.4, 632.0, -610.6, -6.5, -226.8, 12.7);
    Term(2, 0, 1, -4.2, 3.8, -4.1, -4.5, 0.2, 0.6);
    Term(2, -1, 0, 5.3, -0.7, 0.7, 6.1, 0.2, 1.1);
    Term(2, -2, 0, -76.4, -185.1, 260.2, -108.0, 1.6, 0.0);
    Term(2, -3, 0, 66.7, 47.8, -51.4, 69.8, 0.9, 0.3);
    Term(2, -3, 1, 0.6, -1.0, 1.0, 0.6, 0.0, 0.0);
    Term(2, -4, 0, 17.0, 1.4, -1.8, 9.6, 0.0, -0.1);
    Term(2, -5, 0, 1066.2, -518.3, -1.3, -23.9, 1.8, -0.3);
    Term(2, -5, 1, -25.4, -40.3, -0.9, 0.3, 0.0, 0.0);
    Term(2, -5, 2, -0.7, 0.5, 0.0, 0.0, 0.0, 0.0);
    Term(3, 0, 0, -0.1, 28.0, -22.1, -0.2, -12.5, 0.7);
    Term(3, -2, 0, -5.0, -11.5, 11.7, -5.4, 2.1, -1.0);
    Term(3, -3, 0, 16.9, -6.4, 13.4, 26.9, -0.5, 0.8);
    Term(3, -4, 0, 7.2, -13.3, 20.9, 10.5, 0.1, -0.1);
    Term(3, -5, 0, 68.5, 134.3, -166.9, 86.5, 7.1, 15.2);
    Term(3, -5, 1, 3.5, -2.7, 3.4, 4.3, 0.5, -0.4);
    Term(3, -6, 0, 0.6, 1.0, -0.9, 0.5, 0.0, 0.0);
    Term(3, -7, 0, -1.1, 1.7, -0.4, -0.2, 0.0, 0.0);
    Term(4, 0, 0, 0.0, 1.4, -1.0, 0.0, -0.6, 0.0);
    Term(4, -2, 0, -0.3, -0.7, 0.4, -0.2, 0.2, -0.1);
    Term(4, -3, 0, 1.1, -0.6, 0.9, 1.2, 0.1, 0.2);
    Term(4, -4, 0, 3.2, 1.7, -4.1, 5.8, 0.2, 0.1);
    Term(4, -5, 0, 6.7, 8.7, -9.3, 8.7, -1.1, 1.6);
    Term(4, -6, 0, 1.5, -0.3, 0.6, 2.4, 0.0, 0.0);
    Term(4, -7, 0, -1.9, 2.3, -3.2, -2.7, 0.0, -0.1);
    Term(4, -8, 0, 0.4, -1.8, 1.9, 0.5, 0.0, 0.0);
    Term(4, -9, 0, -0.2, -0.5, 0.3, -0.1, 0.0, 0.0);
    Term(4, -10, 0, -8.6, -6.8, -0.4, 0.1, 0.0, 0.0);
    Term(4, -10, 1, -0.5, 0.6, 0.0, 0.0, 0.0, 0.0);
    Term(5, -5, 0, -0.1, 1.5, -2.5, -0.8, -0.1, 0.1);
    Term(5, -6, 0, 0.1, 0.8, -1.6, 0.1, 0.0, 0.0);
    Term(5, -9, 0, -0.5, -0.1, 0.1, -0.8, 0.0, 0.0);
    Term(5, -10, 0, 2.5, -2.2, 2.8, 3.1, 0.1, -0.2);
  end;

  (*sub*)procedure PertUranus; // perturbations by Uranus
  begin
    C[-1] := Cos(M7);
    S[-1] := -Sin(M7);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(1, -1, 0, 0.4, 0.9, 0.0, 0.0, 0.0, 0.0);
    Term(1, -2, 0, 0.4, 0.4, -0.4, 0.3, 0.0, 0.0);
  end;

  (*sub*)procedure PertSUR; // perturbations by Saturn and Uranus
  var
    PHI, X, Y: Double;
  begin
    PHI := (2 * M5 - 6 * M6 + 3 * M7);
    X := Cos(PHI);
    Y := Sin(PHI);
    DL := DL - 0.8 * X + 8.5 * Y;
    DR := DR - 0.1 * X;
    AddThe(X, Y, C5[1], S5[1], X, Y);
    DL := DL + 0.4 * X + 0.5 * Y;
    DR := DR - 0.7 * X + 0.5 * Y;
    DB := DB - 0.1 * X;
  end;

begin // JupiterPos

  DL := 0.0;
  DR := 0.0;
  DB := 0.0;
  M5 := P2 * Frac(0.0565314 + 8.4302963 * T);
  M6 := P2 * Frac(0.8829867 + 3.3947688 * T);
  M7 := P2 * Frac(0.3969537 + 1.1902586 * T);
  C5[0] := 1.0;
  S5[0] := 0.0;
  C5[1] := Cos(M5);
  S5[1] := Sin(M5);
  C5[-1] := C5[1];
  S5[-1] := -S5[1];
  for I := 2 to 5 do
    AddThe(C5[I - 1], S5[I - 1], C5[1], S5[1], C5[I], S5[I]);
  PertSaturn;
  PertUranus;
  PertSUR;
  L := 360.0 * Frac(0.0388910 + M5 / P2 + ((5025.2 + 0.8 * T) * T + DL) / 1296.0E3);
  R := 5.208873 + 0.000041 * T + DR * 1.0E-5;
  B := (227.3 - 0.3 * T + DB) / 3600.0;

end; // JupiterPos

//-----------------------------------------------------------------------

procedure SaturnPos(T: Double; var L, B, R: Double);
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

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I6, I, IT: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C6[I6], S6[I6], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertJupiter; // Kepler terms and perturbations by Jupiter
  var
    I: Integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[1] := Cos(M5);
    S[1] := Sin(M5);
    for I := 0 downto -5 do
      AddThe(C[I], S[I], C[1], -S[1], C[I - 1], S[I - 1]);
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

  (*sub*)procedure PertUranus; // perturbations by Uranus
  var
    I: Integer;
  begin
    C[-1] := Cos(M7);
    S[-1] := -Sin(M7);
    for I := -1 downto - 4 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
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

  (*sub*)procedure PertNep; // perturbations by Neptune
  begin
    C[-1] := Cos(M8);
    S[-1] := -Sin(M8);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(1, -1, 0, -1.3, -1.2, 2.3, -2.5, 0.0, 0.0);
    Term(1, -2, 0, 1.0, -0.1, 0.1, 1.4, 0.0, 0.0);
    Term(2, -2, 0, 1.1, -0.1, 0.2, 3.3, 0.0, 0.0);
  end;

  (*sub*)procedure PertJUR; // perturbations by Jupiter and Uranus
  var
    PHI, X, Y: Double;
  begin
    PHI := (-2 * M5 + 5 * M6 - 3 * M7);
    X := Cos(PHI);
    Y := Sin(PHI);
    DL := DL - 0.8 * X - 0.1 * Y;
    DR := DR - 0.2 * X + 1.8 * Y;
    DB := DB + 0.3 * X + 0.5 * Y;
    AddThe(X, Y, C6[1], S6[1], X, Y);
    DL := DL + (+2.4 - 0.7 * T) * X + (27.8 - 0.4 * T) * Y;
    DR := DR + 2.1 * X - 0.2 * Y;
    AddThe(X, Y, C6[1], S6[1], X, Y);
    DL := DL + 0.1 * X + 1.6 * Y;
    DR := DR - 3.6 * X + 0.3 * Y;
    DB := DB - 0.2 * X + 0.6 * Y;
  end;

begin // SaturnPos

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
    AddThe(C6[I - 1], S6[I - 1], C6[1], S6[1], C6[I], S6[I]);
  PertJupiter;
  PertUranus;
  PertNep;
  PertJUR;
  L := 360.0 * Frac(0.2561136 + M6 / P2 + ((5018.6 + T * 1.9) * T + DL) / 1296.0E3);
  R := 9.557584 - 0.000186 * T + DR * 1.0E-5;
  B := (175.1 - 10.2 * T + DB) / 3600.0;

end; // SaturnPos

//-----------------------------------------------------------------------

procedure UranusPos(T: Double; var L, B, R: Double);
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

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I7, I, IT: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C7[I7], S7[I7], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertJupiter; // perturbations by Jupiter
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(-1, -1, 0, 0.0, 0.0, -0.1, 1.7, -0.1, 0.0);
    Term(0, -1, 0, 0.5, -1.2, 18.9, 9.1, -0.9, 0.1);
    Term(1, -1, 0, -21.2, 48.7, -455.5, -198.8, 0.0, 0.0);
    Term(1, -2, 0, -0.5, 1.2, -10.9, -4.8, 0.0, 0.0);
    Term(2, -1, 0, -1.3, 3.2, -23.2, -11.1, 0.3, 0.1);
    Term(2, -2, 0, -0.2, 0.2, 1.1, 1.5, 0.0, 0.0);
    Term(3, -1, 0, 0.0, 0.2, -1.8, 0.4, 0.0, 0.0);
  end;

  (*sub*)procedure PertSaturn; // perturbations by Saturn
  var
    I: Integer;
  begin
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    for I := -1 downto -3 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
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

  (*sub*)procedure PertNep; // Kepler terms and perturbations by Neptune
  var
    I: Integer;
  begin
    C[-1] := Cos(M8);
    S[-1] := -Sin(M8);
    for I := -1 downto -7 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
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

  (*sub*)procedure PertJSU; // perturbations by Jupiter and Saturn
  var
    I: Integer;
  begin
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    C[-4] := Cos(-4 * M6 + 2 * M5);
    S[-4] := Sin(-4 * M6 + 2 * M5);
    for I := -4 downto -5 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
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

begin // UranusPos

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
    AddThe(C7[I - 1], S7[I - 1], C7[1], S7[1], C7[I], S7[I]);
  for I := 1 to 2 do
  begin
    C7[-I] := C7[I];
    S7[-I] := -S7[I];
  end;
  PertJupiter;
  PertSaturn;
  PertNep;
  PertJSU;
  L := 360.0 * Frac(0.4734843 + M7 / P2 + ((5082.3 + 34.2 * T) * T + DL) / 1296.0E3);
  R := 19.211991 + (-0.000333 - 0.000005 * T) * T + DR * 1.0E-5;
  B := (-130.61 + (-0.54 + 0.04 * T) * T + DB) / 3600.0;

end; // UranusPos

//-----------------------------------------------------------------------

procedure NeptunePos(T: Double; var L, B, R: Double);
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

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I1, I, IT: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  begin
    if IT = 0 then
      AddThe(C8[I1], S8[I1], C[I], S[I], U, V)
    else
    begin
      U := U * T;
      V := V * T
    end;
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertJupiter; // perturbations by Jupiter
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M5);
    S[-1] := -Sin(M5);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(0, -1, 0, 0.1, 0.1, -3.0, 1.8, -0.3, -0.3);
    Term(1, 0, 0, 0.0, 0.0, -15.9, 9.0, 0.0, 0.0);
    Term(1, -1, 0, -17.6, -29.3, 416.1, -250.0, 0.0, 0.0);
    Term(1, -2, 0, -0.4, -0.7, 10.4, -6.2, 0.0, 0.0);
    Term(2, -1, 0, -0.2, -0.4, 2.4, -1.4, 0.4, -0.3);
  end;

  (*sub*)procedure PertSaturn; // perturbations by Saturn
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M6);
    S[-1] := -Sin(M6);
    AddThe(C[-1], S[-1], C[-1], S[-1], C[-2], S[-2]);
    Term(0, -1, 0, -0.1, 0.0, 0.2, -1.8, -0.1, -0.5);
    Term(1, 0, 0, 0.0, 0.0, -8.3, -10.4, 0.0, 0.0);
    Term(1, -1, 0, 13.6, -12.7, 187.5, 201.1, 0.0, 0.0);
    Term(1, -2, 0, 0.4, -0.4, 4.5, 4.5, 0.0, 0.0);
    Term(2, -1, 0, 0.4, -0.1, 1.7, -3.2, 0.2, 0.2);
    Term(2, -2, 0, -0.1, 0.0, -0.2, 2.7, 0.0, 0.0);
  end;

  (*sub*)procedure PertUranus; // Kepler terms and perturbations by Uranus
  var
    I: Integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[-1] := Cos(M7);
    S[-1] := -Sin(M7);
    for I := -1 downto -5 do
      AddThe(C[I], S[I], C[-1], S[-1], C[I - 1], S[I - 1]);
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

begin // NeptunePos

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
    AddThe(C8[I - 1], S8[I - 1], C8[1], S8[1], C8[I], S8[I]);
  PertJupiter;
  PertSaturn;
  PertUranus;
  L := 360.0 * Frac(0.1254046 + M8 / P2 + ((4982.8 - 21.3 * T) * T + DL) / 1296.0E3);
  R := 30.072984 + (0.001234 + 0.000003 * T) * T + DR * 1.0E-5;
  B := (54.77 + (0.26 + 0.06 * T) * T + DB) / 3600.0;

end; // NeptunePos

//-----------------------------------------------------------------------

procedure PlutoPos(T: Double; var L, B, R: Double);

const
  P2 = 6.283185307;
var
  C9, S9: array [0 .. 6] of Double;
  C, S: array [-3 .. 2] of Double;
  M5, M6, M9: Double;
  DL, DR, DB: Double;
  I: Integer;

  (*sub*)function Frac(X: Double): Double;
  begin
    X := X - Trunc(X);
    if (X < 0) then
      X := X + 1.0;
    Frac := X
  end;

  (*sub*)procedure AddThe(C1, S1, C2, S2: Double; var C, S: Double);
  begin
    C := C1 * C2 - S1 * S2;
    S := S1 * C2 + C1 * S2;
  end;

  (*sub*)procedure Term(I9, I: Integer; DLC, DLS, DRC, DRS, DBC, DBS: Double);
  var
    U, V: Double;
  begin
    AddThe(C9[I9], S9[I9], C[I], S[I], U, V);
    DL := DL + DLC * U + DLS * V;
    DR := DR + DRC * U + DRS * V;
    DB := DB + DBC * U + DBS * V;
  end;

  (*sub*)procedure PertJupiter; // Kepler terms and perturbations by Jupiter
  var
    I: Integer;
  begin
    C[0] := 1.0;
    S[0] := 0.0;
    C[1] := Cos(M5);
    S[1] := Sin(M5);
    for I := 0 downto -1 do
      AddThe(C[I], S[I], C[1], -S[1], C[I - 1], S[I - 1]);
    AddThe(C[1], S[1], C[1], S[1], C[2], S[2]);
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

  (*sub*)procedure PertSaturn; // perturbations by Saturn
  var
    I: Integer;
  begin
    C[1] := Cos(M6);
    S[1] := Sin(M6);
    for I := 0 downto -1 do
      AddThe(C[I], S[I], C[1], -S[1], C[I - 1], S[I - 1]);
    Term(1, -1, -29.47, 75.97, -106.4, -204.9, -40.71, -17.55);
    Term(0, 1, -13.88, 18.20, 42.6, -46.1, 1.13, 0.43);
    Term(1, 1, 5.81, -23.48, 15.0, -6.8, -7.48, 3.07);
    Term(2, 1, -10.27, 14.16, -7.9, 0.4, 2.43, -0.09);
    Term(3, 1, 6.86, -10.66, 7.3, -0.3, -2.25, 0.69);
    Term(2, -2, 4.32, 2.00, 0.0, -2.2, -0.24, 0.12);
    Term(1, -2, -5.04, -0.83, -9.2, -3.1, 0.79, -0.24);
    Term(0, -2, 4.25, 2.48, -5.9, -3.3, 0.58, 0.02);
  end;

  (*sub*)procedure PertJUS; // perturbations by Jupiter and Saturn
  var
    PHI, X, Y: Double;
  begin
    PHI := (M5 - M6);
    X := Cos(PHI);
    Y := Sin(PHI);
    DL := DL - 9.11 * X + 0.12 * Y;
    DR := DR - 3.4 * X - 3.3 * Y;
    DB := DB + 0.81 * X + 0.78 * Y;
    AddThe(X, Y, C9[1], S9[1], X, Y);
    DL := DL + 5.92 * X + 0.25 * Y;
    DR := DR + 2.3 * X - 3.8 * Y;
    DB := DB - 0.67 * X - 0.51 * Y;
  end;

  (*sub*)procedure Prec(T: Double; var L, B: Double); // precess. 1950->equinox of date
  const
    DEG = 57.2957795;
  var
    D, PPI, Pis, P, C1, S1, C2, S2, C3, S3, X, Y, Z: Double;
  begin
    D := T + 0.5;
    L := L / DEG;
    B := B / DEG;
    PPI := 3.044;
    Pis := 2.28E-4 * D;
    P := (0.0243764 + 5.39E-6 * D) * D;
    C1 := Cos(Pi);
    C2 := Cos(B);
    C3 := Cos(PPI - L);
    S1 := Sin(Pis);
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

begin // PlutoPos

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
    AddThe(C9[I - 1], S9[I - 1], C9[1], S9[1], C9[I], S9[I]);
  PertJupiter;
  PertSaturn;
  PertJUS;
  L := 360.0 * Frac(0.6232469 + M9 / P2 + DL / 1296.0E3);
  R := 40.7247248 + DR * 1.0E-5;
  B := -3.909434 + DB / 3600.0;
  Prec(T, L, B);

end; // PlutoPos


end.
