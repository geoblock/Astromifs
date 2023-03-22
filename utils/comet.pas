(*-----------------------------------------------------------------------*)
(*                                COMET                                  *)
(*    calculation of unperturbed ephemeris with arbitrary eccentricity   *)
(*                   for comets and minor planets                        *)
(*-----------------------------------------------------------------------*)

program COMET(INPUT,OUTPUT,COMINP);

  uses 
    Astro.Matlib, Astro.Pnulib, Astro.Sphlib, Astro.Sunlib, Astro.Keplib, Astro.Timlib;


  var DAY,MONTH,YEAR,NLINE               : integer;
      D,HOUR,T0,Q,ECC,TEQX0,FAC          : Double;
      MODJD,T,T1,DT,T2,TEQX              : Double;
      X,Y,Z,VX,VY,VZ,XS,YS,ZS            : Double;
      L,B,R,LS,BS,RS,RA,DEC,DELTA,DELTA0 : Double;
      PQR,A,AS                           : Double33;
      COMINP                             : TEXT;


(*-----------------------------------------------------------------------*)
(* GETELM: reads orbital elements from file COMINP                       *)
(*-----------------------------------------------------------------------*)
procedure GETELM (var T0,Q,ECC: double;var PQR:Double33;var TEQX0: double);

  var INC,LAN,AOP: Double;

  begin

    writeln;
    writeln (' COMET: ephemeris calculation for comets and minor planets');
    writeln ('                     version 93/07/01                     ');
    writeln ('        (C) 1993 Thomas Pfleger, Oliver Montenbruck       ');
    writeln;
    writeln (' Orbital elements from file COMINP: ');
    writeln;

    (* open file for reading *)
    (* RESET(COMINP); *)                              (* Standard Pascal *)
    ASSIGN(COMINP,'COMINP.DAT'); RESET(COMINP);       (* TURBO Pascal    *)

    (* display orbital elements *)
    READLN (COMINP,YEAR,MONTH,D);
    writeln('  perihelion time (y m d) ',YEAR:7,MONTH:3,D:6:2);
    READLN(COMINP,Q);
    writeln('  perihelion distance (q) ',  Q:14:7,' AU  ');
    READLN(COMINP,ECC);
    writeln('  eccentricity (e)        ',ECC:14:7);
    READLN(COMINP,INC);
    writeln('  inclination (i)         ',INC:12:5,' deg');
    READLN(COMINP,LAN);
    writeln('  long. of ascending node ',LAN:12:5,' deg');
    READLN(COMINP,AOP);
    writeln('  argument of perihelion  ',AOP:12:5,' deg');
    READLN(COMINP,TEQX0);
    writeln('  equinox                 ',TEQX0:9:2);

    DAY:=Trunc(D); HOUR:=24.0*(D-DAY);
    T0 := ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5) / 36525.0;
    TEQX0 := (TEQX0-2000.0)/100.0;
    GAUSVEC(LAN,INC,AOP,PQR);

  end;

(*-----------------------------------------------------------------------*)
(* GETEPH: reads desired period of time and equinox of the ephemeris     *)
(*-----------------------------------------------------------------------*)
procedure GETEPH(var T1,DT,T2,TEQX: double);

  var YEAR,MONTH,DAY :  Integer;
      EQX,HOUR,JD :     Double;

  begin

    writeln;
    writeln(' Begin and end of the ephemeris: ');
    writeln;
    WRITE  ('  first date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T1 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  final date (yyyy mm dd hh.hhh)            ... ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T2 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    WRITE  ('  step size (dd hh.hh)                      ... ');
    READLN (DAY,HOUR);
    DT :=  ( DAY + HOUR/24.0 ) / 36525.0;
    writeln;
    WRITE  (' Desired equinox of the ephemeris (yyyy.y)  ... ');
    READLN (EQX);
    TEQX := (EQX-2000.0)/100.0;

    writeln; writeln;
    writeln ('    Date      ET   Sun     l      b     r',
             '        Ra          Dec      Distance ');
    writeln (' ':45,'   h  m  s      o  ''  "     (AU) ');

  end;

(*-----------------------------------------------------------------------*)
(* WRTLBR: write L and B in deg,min,sec and R                            *)
(*-----------------------------------------------------------------------*)
procedure WRTLBR(L,B,R: double);
  var H,M: integer;
      S  : Double;
  begin
    DMS(L,H,M,S); WRITE  (H:5,M:3,S:5:1);
    DMS(B,H,M,S); writeln(H:5,M:3,Trunc(S+0.5):3,R:11:6);
  end;
(*------------------------------------------------------------------------*)


begin (* COMET *)


  GETELM (T0,Q,ECC,PQR,TEQX0);  (* read orbital elements              *)
  GETEPH (T1,DT,T2,TEQX);       (* read period of time and equinox    *)

  NLINE := 0;


  PMATECL (TEQX0,TEQX,A);       (* calculate precession matrix         *)

  T := T1;

  repeat

    (* date *)

    MODJD := T*36525.0+51544.5;  CALDAT (MODJD,DAY,MONTH,YEAR,HOUR);

    (* ecliptic coordinates of the sun, equinox TEQX        *)

    Sun200 (T,LS,BS,RS);  CART (RS,BS,LS,XS,YS,ZS);
    PMATECL (T,TEQX,AS);  PRECART (AS,XS,YS,ZS);

    (* heliocentric ecliptic coordinates of the comet       *)

    KEPLER (T0,T,Q,ECC,PQR,X,Y,Z,VX,VY,VZ);
    PRECART (A,X,Y,Z);  PRECART (A,VX,VY,VZ);  POLAR (X,Y,Z,R,B,L);

    (* geometric geocentric coordinates of the comet        *)

    X:=X+XS; Y:=Y+YS; Z:=Z+ZS;  DELTA0 := SQRT ( X*X + Y*Y + Z*Z );

    (* first order correction for light travel time         *)

    FAC:=0.00578*DELTA0;  X:=X-FAC*VX;  Y:=Y-FAC*VY;  Z:=Z-FAC*VZ;
    ECLEQU (TEQX,X,Y,Z);  POLAR (X,Y,Z,DELTA,DEC,RA); RA:=RA/15.0;

    (* output *)

    WRITE(YEAR:4,'/',MONTH:2,'/',DAY:2,HOUR:6:1);
    WRITE(LS:7:1,L:7:1,B:6:1,R:7:3);  WRTLBR(RA,DEC,DELTA0);
    NLINE := NLINE+1; if (NLINE MOD 5) = 0 then writeln;

    (* next time step                                       *)

    T := T + DT;

  until (T2<T);


end.
