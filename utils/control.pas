(*-----------------------------------------------------------------------*)
(*                                CONTROL                                *)
(*              confirmation of orbit determination results              *)
(*-----------------------------------------------------------------------*)

program CONTROL(INPUT,OUTPUT,ORBINP,ORBOUT);

  uses 
    Astro.Matlib, Astro.Pnulib, Astro.Sphlib, Astro.Sunlib, Astro.Keplib, Astro.Timlib;

  var DAY,MONTH,YEAR,I                       : integer;
      D,MODJD,HOUR,T,T0,Q,ECC,TEQX0,TEQX,FAC : Double;
      X,Y,Z,VX,VY,VZ,XS,YS,ZS                : Double;
      L,B,R,LS,BS,RS,RA,DEC,DELTA,DELTA0     : Double;
      PQR,A,AS                               : Double33;
      DATE                                   : array[1..3] of Double;
      ORBINP,ORBOUT                          : TEXT;


(*-----------------------------------------------------------------------*)
(* GETELMF: input orbital elements from file ORBOUT                      *)
(*-----------------------------------------------------------------------*)
procedure GETELMF (var T0,Q,ECC: double;var PQR:Double33;var TEQX0: double);

  var INC,LAN,AOP: Double;
      I          : integer;

  begin

    (* RESET(ORBOUT); *)                              (* Standard Pascal *)
    ASSIGN(ORBOUT,'ORBOUT.DAT'); RESET(ORBOUT);       (* Turbo Pascal    *)

    READLN(ORBOUT,YEAR,MONTH,D); DAY:=Trunc(D); HOUR:=24.0*(D-DAY);
    writeln(' Perihelion date (y m d) ',YEAR:4,'/',MONTH:2,'/',D:6:2);
    T0 := ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5) / 36525.0;
    READLN(ORBOUT,Q);    READLN(ORBOUT,ECC);  READLN(ORBOUT,INC);
    READLN(ORBOUT,LAN);  READLN(ORBOUT,AOP);  READLN(ORBOUT,TEQX0);
    writeln(' Perihelion distance (q) ',Q:12:7,' AU');
    writeln(' Eccentricity (e)        ',ECC:12:7);
    writeln(' Inclination (i)         ',INC:10:5,' deg');
    writeln(' Long. of ascending node ',LAN:10:5,' deg');
    writeln(' Argument of perihelion  ',AOP:10:5,' deg');
    writeln(' Equinox                 ',TEQX0:7:2);

    TEQX0 := (TEQX0-2000)/100;
    GAUSVEC(LAN,INC,AOP,PQR);

  end;

(*-----------------------------------------------------------------------*)
procedure WRTLBR(L,B,R: double);
  var H,M: integer;
      S  : Double;
  begin
    DMS(L,H,M,S); WRITE  (H:5,M:3,S:5:1);
    DMS(B,H,M,S); writeln(H:5,M:3,Trunc(S+0.5):3,R:11:6);
  end;
(*-----------------------------------------------------------------------*)

begin

  (* input                                                  *)

  GETELMF (T0,Q,ECC,PQR,TEQX0); (* orbital elements         *)

  (* RESET(ORBINP); *)                               (* Standard Pascal *)
  ASSIGN(ORBINP,'ORBINP.DAT'); RESET(ORBINP);        (* Turbo Pascal    *)


  READLN(ORBINP);
  for I := 1 to 3 do
    begin
      READLN(ORBINP,YEAR,MONTH,DAY,HOUR);
      DATE[I] := MJD(DAY,MONTH,YEAR,HOUR);
    end;
  READ(ORBINP,TEQX); TEQX:=(TEQX-2000.0)/100.0;

  writeln;
  writeln ('    Date      ET   Sun     l      b     r',
           '        Ra          Dec      Distance ');
  writeln (' ':45,'   h  m  s      o  ''  "    (AU) ');

  PMATECL (TEQX0,TEQX,A);

  for I:=1 to 3 do

   begin

    (* date                                                  *)

    MODJD := DATE[I];
    CALDAT (MODJD,DAY,MONTH,YEAR,HOUR);
    T := (MODJD-51544.5)/36525.0;

    (* ecliptic coordinates of the sun, equinox TEQX         *)

    Sun200 (T,LS,BS,RS);
    CART (RS,BS,LS,XS,YS,ZS);
    PMATECL (T,TEQX,AS); PRECART (AS,XS,YS,ZS);

    (* heliocentric ecliptic coordinates of the comet/planet *)

    KEPLER (T0,T,Q,ECC,PQR,X,Y,Z,VX,VY,VZ);
    PRECART (A,X,Y,Z);   PRECART (A,VX,VY,VZ);
    POLAR (X,Y,Z,R,B,L);

    (* geometric geocentric coordinates of the comet         *)

    X:=X+XS; Y:=Y+YS; Z:=Z+ZS;
    DELTA0 := SQRT ( X*X + Y*Y + Z*Z );

    (* first order correction for light travel time          *)

    FAC := 0.00578*DELTA0;
    X := X-FAC*VX;  Y := Y-FAC*VY;  Z := Z-FAC*VZ;
    ECLEQU (TEQX,X,Y,Z);
    POLAR (X,Y,Z,DELTA,DEC,RA); RA:=RA/15.0;

    (* output                                               *)

    WRITE(YEAR:5,'/',MONTH:2,'/',DAY:2,HOUR:5:1);
    WRITE(LS:7:1,L:7:1,B:6:1,R:7:3);
    WRTLBR(RA,DEC,DELTA0);

   end;


end.
