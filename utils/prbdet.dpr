(*-----------------------------------------------------------------------*)
(*                               ORBDET                                  *)
(*         Gaussian orbit determination from three observations          *)
(*               using the abbreviated method of Bucerius                *)
(*-----------------------------------------------------------------------*)

program ORBDET(Input,Output,ORBINP,ORBOUT);

{$APPTYPE CONSOLE}

  uses 
    Astro.Matlib, Astro.Pnulib, Astro.Sphlib, 
    Astro.Sunlib, Astro.Timlib, Astro.Keplib;

  type CHAR80 = array[1..80] of CHAR;

  var  TEQX                : Double;
       TP,Q,ECC,INC,LAN,AOP: Double;
       JD0                 : REAL3;
       RSUN,E              : MAT3X;
       HEADER              : CHAR80;
       ORBINP,ORBOUT       : TEXT;


(*-----------------------------------------------------------------------*)
(* START: reads the input file and preprocesses the observational data   *)
(*                                                                       *)
(* output:                                                               *)
(*   RSUN:  matrix of three Sun position vectors in ecliptic coordinates *)
(*   E:     matrix of three observation direction unit vectors           *)
(*   JD:    julian date of the three observation times                   *)
(*   TEQX:  equinox of RSUN and E (in Julian centuries since J2000)      *)
(*-----------------------------------------------------------------------*)

procedure START (var HEADER: CHAR80;
                 var RSUN,E: MAT3X; var JD0: REAL3; var TEQX: Double);

  var DAY,MONTH,YEAR,D,M,I    :  Integer;
      UT,S,DUMMY              :  Double;
      EQX0,EQX,TEQX0          :  Double;
      LS,BS,RS,LP,BP,RA,DEC,T :  REAL3;
      A,ASI                    :  Double33;
      ORBINP                  :  TEXT;

  begin

    (* open input file                                                   *)

    (* RESET(ORBINP); *)                              (* Standard Pascal *)
    ASSIGN(ORBINP,'ORBINP.DAT'); RESET(ORBINP);       (* Turbo Pascal    *)

    (* read data from file ORBINP                                        *)

    for I:=1 to 80 do                                 (* header          *)
      if NOT(EOLN(ORBINP)) then read(ORBINP,HEADER[I]) else HEADER[I]:=' ';
    readln(ORBINP);
    for I := 1 to 3 do                                (* 3 observations  *)
      begin
        READ  (ORBINP,YEAR,MONTH,DAY,UT);                        (* date *)
        READ  (ORBINP,D,M,S); DDD(D,M,S,RA[I]);                  (* RA   *)
        readln(ORBINP,D,M,S); DDD(D,M,S,DEC[I]);                 (* Dec  *)
        RA[I]:=RA[I]*15.0;
        JD0[I] := 2400000.5+MJD(DAY,MONTH,YEAR,UT);
        T[I]   := (JD0[I]-2451545.0)/36525.0;
      end;
    writeln;
    readln(ORBINP,EQX0); TEQX0:=(EQX0-2000.0)/100.0;  (* equinox         *)

    (* desired equinox of the orbital elements                           *)

    read(ORBINP,EQX ); TEQX :=(EQX -2000.0)/100.0;

    (* calculate initial data of the orbit determination                 *)

    PMATECL(TEQX0,TEQX,A);
    for I := 1 to 3 do
      begin
        CART   (1.0,DEC[I],RA[I],E[I,X],E[I,Y],E[I,Z]);
        EQUECL (     TEQX0      ,E[I,X],E[I,Y],E[I,Z]);
        PRECART(      A         ,E[I,X],E[I,Y],E[I,Z]);
        POLAR  (E[I,X],E[I,Y],E[I,Z],DUMMY,BP[I],LP[I]);
        PMATECL( T[I],TEQX,ASI);
        Sun200 (T[I],LS[I],BS[I],RS[I]);
        CART   (RS[I],BS[I],LS[I],RSUN[I,X],RSUN[I,Y],RSUN[I,Z]);
        PRECART(      ASI         ,RSUN[I,X],RSUN[I,Y],RSUN[I,Z]);
      end;

    writeln('   ORBDET: orbit determination from three observations ');
    writeln('                     version 93/07/01                  ');
    writeln('        (c) 1993 Thomas Pfleger, Oliver Montenbruck    ');
    writeln; writeln;
    writeln(' Summary of orbit determination ');
    writeln;
    write  ('  '); for I:=1 to 78 do write(HEADER[I]); writeln;
    writeln;
    writeln(' Initial data (ecliptic geocentric coordinates (in deg))');
    writeln;
    writeln('  Julian Date        ', JD0[1]:12:2,JD0[2]:12:2,JD0[3]:12:2);
    writeln('  Solar longitude    ',  LS[1]:12:2, LS[2]:12:2, LS[3]:12:2);
    writeln('  Planet/Comet Longitude',LP[1]:9:2, LP[2]:12:2, LP[3]:12:2);
    writeln('  Planet/Comet Latitude ',BP[1]:9:2, BP[2]:12:2, BP[3]:12:2);
    writeln; writeln;

  end;

(*-----------------------------------------------------------------------*)
(* DUMPELEM: output of orbital elements (screen)                         *)
(*-----------------------------------------------------------------------*)
procedure DUMPELEM(TP,Q,ECC,INC,LAN,AOP,TEQX: double);
  var DAY,MONTH,YEAR: integer;
      MODJD,UT      : Double;
  begin
    MODJD := TP*36525.0 + 51544.5;
    CALDAT( MODJD, DAY,MONTH,YEAR,UT);
    writeln(' Orbital elements',
            ' (Equinox ','J',100.0*TEQX+2000.0:8:2,')');
    writeln;
    writeln('  Perihelion date      tp    ',
            YEAR:4,'/',MONTH:2,'/',DAY:2,UT:8:4,'h',
            '  (JD',MODJD+2400000.5:11:2,')');
    writeln('  Perihelion distance  q[AU] ',  Q:12:6);
    writeln('  Semi-major axis      a[AU] ',  Q/(1-ECC):12:6);
    writeln('  Eccentricity         e     ',  ECC:12:6);
    writeln('  Inclination          i     ',  INC:10:4,' degrees');
    writeln('  Ascending node       Omega ',  LAN:10:4,' degrees');
    writeln('  Long. of perihelion  pi    ',  AOP+LAN:10:4,' degrees');
    writeln('  Arg. of perihelion   omega ',  AOP:10:4,' degrees');
    writeln;
  end;

(*-----------------------------------------------------------------------*)
(* SAVEELEM: output of orbital elements (file)                           *)
(*-----------------------------------------------------------------------*)
procedure SAVEELEM(TP,Q,ECC,INC,LAN,AOP,TEQX: double;HEADER: CHAR80);

  var I,DAY,MONTH,YEAR: integer;
      MODJD,UT        : Double;

  begin

    (* open file for writing *)

    (* REWRITE(ORBOUT); *)                            (* Standard Pascal *)
    ASSIGN(ORBOUT,'ORBOUT.DAT'); REWRITE(ORBOUT);     (* Turbo Pascal    *)

    MODJD := TP*36525.0 + 51544.5;
    CALDAT( MODJD, DAY,MONTH,YEAR,UT);
    write  (ORBOUT,YEAR:5,MONTH:3,(DAY+UT/24.0):7:3,'!':6);
    writeln(ORBOUT,' perihelion time T0 (y m d.d)  =  JD ',
                   (MODJD+2400000.5):12:3);
    writeln(ORBOUT, Q :12:6,'!': 9,' q  ( a =',Q/(1-ECC):10:6,' )');
    writeln(ORBOUT,ECC:12:6,'!': 9,' e ');
    writeln(ORBOUT,INC:10:4,'!':11,' i ');
    writeln(ORBOUT,LAN:10:4,'!':11,' long.asc.node ');
    writeln(ORBOUT,AOP:10:4,'!':11,
                   ' arg.perih. ( long.per. = ',AOP+LAN:9:4,' )');
    writeln(ORBOUT,TEQX*100.0+2000.0:8:2,'!':13,' equinox (J)');
    write  (ORBOUT,'! ');
    for I:=1 to 78 do write(ORBOUT,HEADER[I]);

    RESET(ORBOUT); (* close file *)

  end;

(*-----------------------------------------------------------------------*)
(* RETARD: light-time correction                                         *)
(*   JD0: times of observation (t1',t2',t3') (Julian Date)               *)
(*   RHO: three geocentric distances (in AU)                             *)
(*   JD:  times of light emittance (t1,t2,t3) (Julian Date)              *)
(*   TAU: scaled time differences                                        *)
(*-----------------------------------------------------------------------*)
procedure RETARD ( JD0,RHO: REAL3; var JD,TAU: REAL3);
  const KGAUSS = 0.01720209895;  A = 0.00578;
  var   I: integer;
  begin
    for I:=1 to 3 do  JD[I]:=JD0[I]-A*RHO[I];
    TAU[1] := KGAUSS*(JD[3]-JD[2]);  TAU[2] := KGAUSS*(JD[3]-JD[1]);
    TAU[3] := KGAUSS*(JD[2]-JD[1]);
  end;

(*-----------------------------------------------------------------------*)
(* GAUSS: iteration of the abbreviated Gauss method                      *)
(*                                                                       *)
(*  RSUN: three vectors of geocentric Sun positions                      *)
(*  E   : three unit vectors of geocentric observation directions        *)
(*  JD0 : three observation times (Julian Date)                          *)
(*  TP  : time of perihelion passage (Julian centuries since J2000)      *)
(*  Q   : perihelion distance                                            *)
(*  ECC : eccentricity                                                   *)
(*  INC : inclination                                                    *)
(*  LAN : longitude of the ascending node                                *)
(*  AOP : argument of perihelion                                         *)
(*-----------------------------------------------------------------------*)

procedure GAUSS ( RSUN,E: MAT3X; JD0:REAL3;
                  var TP,Q,ECC,INC,LAN,AOP: Double);

  const EPS_RHO =1.0E-8;

  var I,J              : integer;
      S                : INDEX;
      RHOOLD,DET       : Double;
      JD,RHO,N,TAU,ETA : REAL3;
      DI               : VECTOR;
      RPL              : MAT3X;
      DD               : Double33;

  begin

    (* calculate initial approximations of n1 and n3 *)

    N[1] := (JD0[3]-JD0[2]) / (JD0[3]-JD0[1]);     N[2] := -1.0;
    N[3] := (JD0[2]-JD0[1]) / (JD0[3]-JD0[1]);

    (* calculate matrix D and its determinant (det(D) = e3.d3) *)

    CROSS(E[2],E[3],DI);  for J:=1 to 3 do DD[1,J]:=DOT(DI,RSUN[J]);
    CROSS(E[3],E[1],DI);  for J:=1 to 3 do DD[2,J]:=DOT(DI,RSUN[J]);
    CROSS(E[1],E[2],DI);  for J:=1 to 3 do DD[3,J]:=DOT(DI,RSUN[J]);
    DET := DOT(E[3],DI);

    writeln; writeln(' Iteration of the geocentric distances rho [AU] ');
    writeln;

    RHO[2] := 0;

    (* Iterate until distance rho[2] does not change any more *)

    RHO[2] := 0;

    repeat

       RHOOLD := RHO[2];

      (* geocentric distance rho from n1 and n3 *)
      for I := 1 to 3 do
        RHO[I]:=( N[1]*DD[I,1] - DD[I,2] + N[3]*DD[I,3] ) / (N[I]*DET);

      (* apply light-time correction and calculate time differences *)
      RETARD (JD0,RHO,JD,TAU);

      (* heliocentric coordinate vectors *)
      for I := 1 to 3 do
        for S := X to Z do
          RPL[I,S] := RHO[I]*E[I,S]-RSUN[I,S];

      (* sector/triangle ratios eta[i] *)
      ETA[1] := FIND_ETA( RPL[2], RPL[3], TAU[1] );
      ETA[2] := FIND_ETA( RPL[1], RPL[3], TAU[2] );
      ETA[3] := FIND_ETA( RPL[1], RPL[2], TAU[3] );

      (* improvement of the sector/triangle ratios *)
      N[1] := ( TAU[1]/ETA[1] ) / (TAU[2]/ETA[2]);
      N[3] := ( TAU[3]/ETA[3] ) / (TAU[2]/ETA[2]);
      writeln('  rho',' ':16,RHO[1]:12:8,RHO[2]:12:8,RHO[3]:12:8);

    until ( abs(RHO[2]-RHOOLD) < EPS_RHO );

    writeln; writeln(' Heliocentric distances [AU]:'); writeln;
    writeln('  r  ',' ':16,
            NORM(RPL[1]):12:8,NORM(RPL[2]):12:8,NORM(RPL[3]):12:8);
    writeln; writeln;

    (* derive orbital elements from first and third observation *)

    ELEMENT ( JD[1],JD[3],RPL[1],RPL[3], TP,Q,ECC,INC,LAN,AOP );


  end;

(*----------------------------------------------------------------------*)

begin

  START(HEADER,RSUN,E,JD0,TEQX);

  GAUSS(RSUN,E,JD0,TP,Q,ECC,INC,LAN,AOP);

  DUMPELEM(TP,Q,ECC,INC,LAN,AOP,TEQX);
  SAVEELEM(TP,Q,ECC,INC,LAN,AOP,TEQX,HEADER);

  (* check solution  *)

  writeln;
  if (DOT(E[2],RSUN[2])>0) then
    writeln (' Warning: observation in hemisphere of conjunction;',
             '  possible second solution');
  if (ECC>1.1) then
    writeln (' Warning: probably not a realistic solution (e>1.1) ');
  if ( (abs(Q-0.985)<0.1) and (abs(ECC-0.015)<0.05) ) then
    writeln (' Warning: probably Earth''s orbit solution');

end.

(*-----------------------------------------------------------------------*)

