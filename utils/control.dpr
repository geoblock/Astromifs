(*-----------------------------------------------------------------------*)
(*                                CONTROL                                *)
(*              confirmation of orbit determination results              *)
(*-----------------------------------------------------------------------*)

program CONTROL(Input,Output,OrbInput,OrbOutput);

{$APPTYPE CONSOLE}

  uses 
    Astro.Matlib,
    Astro.Pnulib,
    Astro.Sphlib,
    Astro.Sunlib,
    Astro.Kepler,
    Astro.Timlib;

  var DAY,MONTH,YEAR,I                       : integer;
      D,MODJD,HOUR,T,T0,Q,ECC,TEQX0,TEQX,FAC : Double;
      X,Y,Z,VX,VY,VZ,XS,YS,ZS                : Double;
      L,B,R,LS,BS,RS,RA,DEC,DELTA,DELTA0     : Double;
      PQR,A,ASI                               : Double33;
      DATE                                   : array[1..3] of Double;
      OrbInput,OrbOutput                          : TEXT;


(*-----------------------------------------------------------------------*)
(* GETELMF: input orbital elements from file OrbOutput                      *)
(*-----------------------------------------------------------------------*)
procedure GETELMF (var T0,Q,ECC: double;var PQR:Double33;var TEQX0: double);

  var INC,LAN,AOP: Double;
      I          : integer;

  begin

    (* RESET(OrbOutput); *)                              (* Standard Pascal *)
    ASSIGN(OrbOutput,'OrbOutput.DAT'); RESET(OrbOutput);       (* Turbo Pascal    *)

    readln(OrbOutput,YEAR,MONTH,D); DAY:=Trunc(D); HOUR:=24.0*(D-DAY);
    writeln(' Perihelion date (y m d) ',YEAR:4,'/',MONTH:2,'/',D:6:2);
    T0 := ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5) / 36525.0;
    readln(OrbOutput,Q);    readln(OrbOutput,ECC);  readln(OrbOutput,INC);
    readln(OrbOutput,LAN);  readln(OrbOutput,AOP);  readln(OrbOutput,TEQX0);
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
    DMS(L,H,M,S); write  (H:5,M:3,S:5:1);
    DMS(B,H,M,S); writeln(H:5,M:3,Trunc(S+0.5):3,R:11:6);
  end;
(*-----------------------------------------------------------------------*)

begin

  (* input                                                  *)

  GETELMF (T0,Q,ECC,PQR,TEQX0); (* orbital elements         *)

  (* RESET(OrbInput); *)                               (* Standard Pascal *)
  ASSIGN(OrbInput,'OrbInput.DAT'); RESET(OrbInput);        (* Turbo Pascal    *)


  readln(OrbInput);
  for I := 1 to 3 do
    begin
      readln(OrbInput,YEAR,MONTH,DAY,HOUR);
      DATE[I] := MJD(DAY,MONTH,YEAR,HOUR);
    end;
  read(OrbInput,TEQX); TEQX:=(TEQX-2000.0)/100.0;

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
    PMATECL (T,TEQX,ASI); PRECART (ASI,XS,YS,ZS);

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

    write(YEAR:5,'/',MONTH:2,'/',DAY:2,HOUR:5:1);
    write(LS:7:1,L:7:1,B:6:1,R:7:3);
    WRTLBR(RA,DEC,DELTA0);

   end;


end.
