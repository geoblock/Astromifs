(*-----------------------------------------------------------------------*)
(*                                                                       *)
(*                                NUMINT                                 *)
(*        Numerical integration of perturbed minor planet orbits         *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

program NUMINT ( INPUT,OUTPUT, NUMINP );

  USES  
    Astro.Matlib,Astro.Timlib,Astro.Sphlib, 
    Astro.Pnulib, Astro.Keplib, Astro.Sunlib, 
    Astro.Delib, Astro.Plalib;

  const J2000  = 0.0;

  var DAY,MONTH,YEAR,NLINE,IFLAG         : integer;
      D,HOUR,T_EPOCH, TEQX0,TEQX, FAC    : Double;
      A,E,M,INC,LAN,AOP                  : Double;
      MJD_START,MJD_END, T,T1,DT,T2      : Double;
      XX,YY,ZZ, VX,VY,VZ, XS,YS,ZS       : Double;
      L,B,R,LS,BS,RS,RA,DEC,DELTA,DELTA0 : Double;
      PQR                                : Double33;
      EQX0_TO_J2000, J2000_TO_EQX, AS    : Double33;
      Y                                  : DE_EQN_VECTOR;      
      WORK                               : DE_WORKSPACE_RECORD; 
      NUMINP                             : TEXT;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* WRTELM: writes orbital elements                                       *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

procedure WRTELM ( YEAR,MONTH: integer; D: Double;
                   A,E,INC,LAN,AOP,M, TEQX: Double);
  begin
    writeln('  Epoch (y m d)           ',YEAR:7,MONTH:3,D:6:2);
    writeln('  Semi-major axis (a)     ',  A:14:7,' AU  ');
    writeln('  Eccentricity (e)        ',  E:14:7);
    writeln('  Inclination (i)         ',INC:12:5,' deg');
    writeln('  Long. of ascending node ',LAN:12:5,' deg');
    writeln('  Argument of perihelion  ',AOP:12:5,' deg');
    writeln('  Mean anomaly (M)        ',  M:12:5,' deg');
    writeln('  Equinox                 ',2000.0+100.0*TEQX:9:2);
  end;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* GETELM: reads orbital elements from file NUMINP                       *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

procedure GETELM (var T,A,E,M: Double; var PQR: Double33; var TEQX0: Double);

  var YEAR,MONTH,DAY : integer;
      D,HOUR         : Double;
      INC,LAN,AOP    : Double;

  begin

    writeln;
    writeln ('                          NUMINT                          ');
    writeln ('  Numerical integration of perturbed minor planet orbits  ');
    writeln ('                     version 93/07/01                     ');
    writeln ('        (c) 1993 Thomas Pfleger, Oliver Montenbruck       ');
    writeln;
    writeln (' Orbital elements from file NUMINP: ');
    writeln;

    (* Open file for reading *)

    (* RESET(NUMINP); *)                         (* Standard Pascal *)
    ASSIGN(NUMINP,'NUMINP.DAT'); RESET(NUMINP);  (* TURBO Pascal *)

    (* Display orbital elements *)

    READLN (NUMINP,YEAR,MONTH,D);
    READLN (NUMINP,A);
    READLN (NUMINP,E);
    READLN (NUMINP,INC);
    READLN (NUMINP,LAN);
    READLN (NUMINP,AOP);
    READLN (NUMINP,M);
    READLN (NUMINP,TEQX0);

    DAY:=Trunc(D); HOUR:=24.0*(D-DAY);
    T := ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5) / 36525.0;
    TEQX0 := (TEQX0-2000.0)/100.0;
    GAUSVEC(LAN,INC,AOP,PQR);

    WRTELM (YEAR,MONTH,D, A,E,INC,LAN,AOP,M, TEQX0 );

  end;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* GETEPH: reads desired period of time and equinox of the ephemeris     *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

procedure GETEPH ( var T1,DT,T2,TEQX: Double);

  var YEAR,MONTH,DAY :  Integer;
      EQX,HOUR,JD :     Double;

  begin

    writeln;
    writeln(' Begin and end of the ephemeris: ');
    writeln;
    WRITE  ('  First date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T1 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  Final date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T2 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  Step size (dd hh.hh)                      ... ');
    READLN (DAY,HOUR);
    DT :=  ( DAY + HOUR/24.0 ) / 36525.0;
    writeln;
    WRITE  (' Desired equinox of the ephemeris (yyyy.y)  ... ');
    READLN (EQX);
    TEQX := (EQX-2000.0)/100.0;

  end;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* WRTLBR: write L and B in deg,min,sec and R                            *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

procedure WRTLBR(L,B,R: double);
  var H,M: integer;
      S  : Double;
  begin
    DMS(L,H,M,S); WRITE  (H:5,M:3,S:5:1);
    DMS(B,H,M,S); writeln(H:5,M:3,Trunc(S+0.5):3,R:11:6);
  end;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* ACCEL: computes the acceleration                                      *)
(*                                                                       *)
(*   T         Time in Julian centuries since J2000                      *)
(*   X,Y,Z     Heliocentric ecliptic coordinates (in AU)                 *)
(*   AX,AY,AZ  Acceleration (in AU/d**2)                                 *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

procedure ACCEL ( T, X,Y,Z: Double; var AX,AY,AZ: Double);

  const K_GAUSS =  0.01720209895;     (* Gaussian gravitational constant *)

  var   Planet                :  PlanetType;
        GM_SUN, R_SQR,R,MU_R3 :  Double;
        XP,YP,ZP, DX,DY,DZ    :  Double;
        MU                    :  array[PlanetType] of Double;

  begin

    (* Grav. constant * solar and planetary masses in AU**3/d**2 *)

    GM_SUN := K_GAUSS*K_GAUSS;

    MU[Mercury] := GM_SUN / 6023600.0;
    MU[Venus]   := GM_SUN /  408523.5;
    MU[Earth]   := GM_SUN /  328900.5;
    MU[Mars]    := GM_SUN / 3098710.0;
    MU[Jupiter] := GM_SUN /    1047.355;
    MU[Saturn]  := GM_SUN /    3498.5;
    MU[Uranus]  := GM_SUN /   22869.0;
    MU[Neptune] := GM_SUN /   19314.0;
    MU[Pluto]   := GM_SUN / 3000000.0;

    (* Solar attraction *)

    R_SQR := ( X*X + Y*Y + Z*Z );  R := Sqrt(R_SQR);
    MU_R3 := GM_SUN / (R_SQR*R);
    AX := -MU_R3*X;  AY := -MU_R3*Y;  AZ := -MU_R3*Z;

    (* Planetary perturbation *)

    for Planet := Mercury to Pluto do

      begin

        (* Planetary coordinates *)

        POSITION ( Planet, T, XP,YP,ZP );
        DX:=X-XP;  DY:=Y-YP;  DZ:=Z-ZP;

        (* Direct acceleration   *)

        R_SQR := ( DX*DX + DY*DY + DZ*DZ );  R := Sqrt(R_SQR);
        MU_R3 := MU[Planet] / (R_SQR*R);
        AX := AX-MU_R3*DX;  AY := AY-MU_R3*DY;  AZ := AZ-MU_R3*DZ;

        (* Indirect acceleration *)

        R_SQR := ( XP*XP + YP*YP + ZP*ZP );  R := Sqrt(R_SQR);
        MU_R3 := MU[Planet] / (R_SQR*R);
        AX := AX-MU_R3*XP;  AY := AY-MU_R3*YP;  AZ := AZ-MU_R3*ZP;

      end;

  end;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* F: computes the time derivative of the state vector                   *)
(*                                                                       *)
(*   X   Time (Modified Julian Date)                                     *)
(*   Y   State vector (x,y,z in AU, vx,vy,vz in AU/d)                    *)
(*   YP  Derivative (vx,vy,vz in AU/d, ax,ay,az in AU/d**2               *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

{$F+ Turbo Pascal compiler directive: force far call}

procedure F ( X: Double; Y: DE_EQN_VECTOR; var DYDX: DE_EQN_VECTOR );

  var T : Double;

  begin

    (* Time in Julian centuries since J2000 *)

    T := ( X-51544.5 ) / 36525.0;

    (* Derivative of the state vector *)

    DYDX[1]:=Y[4]; DYDX[2]:=Y[5]; DYDX[3]:=Y[6];

    ACCEL ( T, Y[1],Y[2],Y[3], DYDX[4],DYDX[5],DYDX[6] );

  end;

{$F- Turbo Pascal compiler directive: end far call}


(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* INTEGRATE: Integrates the equation of motion                          *)
(*                                                                       *)
(*   Y        State vector (x,y,z in AU, vx,vy,vz in AU/d)               *)
(*   MJD      Epoch (Modified Julian Date)                               *)
(*   MJD_END  Final epoch (Modified Julian Date)                         *)
(*   IFLAG    Return code                                                *)
(*   WORK     Work space                                                 *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

procedure INTEGRATE ( var Y            : DE_EQN_VECTOR;
                      var MJD, MJD_END : Double;
                      var IFLAG        : integer;
                      var WORK         : DE_WORKSPACE_RECORD );

  const EPS = 1.0E-8;  (* Accuracy *)

  var   RELERR, ABSERR: Double;

  begin

    RELERR := EPS;    (* Relative accuracy requirement *)
    ABSERR := 0.0;    (* Absolute accuracy requirement *)

    if ( MJD_END <> MJD ) then
      begin
        repeat
          DE ( F, 6, Y, MJD, MJD_END, RELERR,ABSERR, IFLAG, WORK );
        until ((abs(IFLAG)=2) or (IFLAG=6));
        if (IFLAG=6) then writeln (' Illegal input in DE ');
      end;

  end;

(*-----------------------------------------------------------------------*)


begin (* NUMINT *)


  (* Read orbital elements, prediction interval and equinox *)

  GETELM (T_EPOCH,A,E,M,PQR,TEQX0);
  GETEPH (T1,DT,T2,TEQX);

  (* Calculate precession matrices *)

  PMATECL (TEQX0,J2000,EQX0_TO_J2000);
  PMATECL (J2000,TEQX,J2000_TO_EQX);

  (* Heliocentric position and velocity vector at epoch *)
  (* referred to ecliptic and equinox of J2000          *)

  ELLIP  ( M,A,E, XX,YY,VX,VY );
  ORBECL ( XX,YY,PQR, Y[1],Y[2],Y[3] );
  ORBECL ( VX,VY,PQR, Y[4],Y[5],Y[6] );

  PRECART ( EQX0_TO_J2000,Y[1],Y[2],Y[3] );
  PRECART ( EQX0_TO_J2000,Y[4],Y[5],Y[6] );


  (* Start integration: propagate state vector from epoch *)
  (* to start of ephemeris                                *)

  MJD_START := T_EPOCH*36525.0 + 51544.5;
  MJD_END   := T1*36525.0 + 51544.5;

  IFLAG     := 1;   (* Initialization flag *)

  INTEGRATE ( Y, MJD_START, MJD_END, IFLAG, WORK );


  (* Orbital elements at start of ephemeris *)

  XX:=Y[1]; YY:=Y[2]; ZZ:=Y[3];    (* Copy J2000 state vector *)
  VX:=Y[4]; VY:=Y[5]; VZ:=Y[6];

  PRECART (J2000_TO_EQX,XX,YY,ZZ);              (* Precession *)
  PRECART (J2000_TO_EQX,VX,VY,VZ);

  XYZKEP ( XX,YY,ZZ,VX,VY,VZ, A,E,INC,LAN,AOP,M ); (* Convert *)

  writeln; writeln;
  writeln (' Orbital elements at start epoch:');
  writeln;

  CALDAT (MJD_END,DAY,MONTH,YEAR,HOUR);

  WRTELM ( YEAR,MONTH,DAY+HOUR/24.0,        (* Print elements *)
           A,E,INC,LAN,AOP,M, TEQX );


  (* Create ephemeris *)

  writeln; writeln;
  writeln ('    Date      ET   Sun     l      b     r',
           '        RA          Dec      Distance ');
  writeln (' ':45,'   h  m  s      o  ''  "     (AU) ');

  NLINE := 0;
  T     := T1;

  repeat

    (* Integrate orbit to time T *)

    MJD_END := T*36525.0 + 51544.5;

    INTEGRATE ( Y, MJD_START, MJD_END, IFLAG, WORK );

    (* Heliocentric ecliptic coordinates, equinox TEQX *)

    XX:=Y[1]; YY:=Y[2]; ZZ:=Y[3];    (* Copy J2000 state vector *)
    VX:=Y[4]; VY:=Y[5]; VZ:=Y[6];

    PRECART (J2000_TO_EQX,XX,YY,ZZ);              (* Precession *)
    PRECART (J2000_TO_EQX,VX,VY,VZ);

    POLAR (XX,YY,ZZ,R,B,L);


    (* Ecliptic coordinates of the Sun, equinox TEXQ *)

    Sun200 (T,LS,BS,RS);  CART (RS,BS,LS,XS,YS,ZS);
    PMATECL (T,TEQX,AS);  PRECART (AS,XS,YS,ZS);


    (* Geometric geocentric coordinates *)

    XX:=XX+XS; YY:=YY+YS; ZZ:=ZZ+ZS;


    (* First order correction for light travel time *)

    DELTA0 := SQRT ( XX*XX + YY*YY + ZZ*ZZ );
    FAC    := 0.00578*DELTA0;

    XX:=XX-FAC*VX;  YY:=YY-FAC*VY;  ZZ:=ZZ-FAC*VZ;

    ECLEQU (TEQX,XX,YY,ZZ);
    POLAR (XX,YY,ZZ,DELTA,DEC,RA); RA:=RA/15.0;


    (* Output *)

    CALDAT (MJD_END,DAY,MONTH,YEAR,HOUR);

    WRITE(YEAR:4,'/',MONTH:2,'/',DAY:2,HOUR:6:1);
    WRITE(LS:7:1,L:7:1,B:6:1,R:7:3);  WRTLBR(RA,DEC,DELTA0);
    NLINE := NLINE+1; if (NLINE MOD 5) = 0 then writeln;


    (* Next time step *)

    T := T + DT;


  until (T2<T);


end.


