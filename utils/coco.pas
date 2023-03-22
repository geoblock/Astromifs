(*-----------------------------------------------------------------------*)
(*                                  COCO                                 *)
(*                         coordinate conversion                         *)
(*-----------------------------------------------------------------------*)

program COCO(INPUT,OUTPUT);

  uses 
    Astro.Matlib, Astro.Pnulib, Astro.Sphlib, Astro.Sunlib, Astro.Timlib;

  var X,Y,Z,XS,YS,ZS: Double;
      T,TEQX,TEQXN  : Double;
      LS,BS,RS      : Double;
      A             : Double33;
      ECLIPT        : Boolean;
      MODE          : CHAR;


(*-----------------------------------------------------------------------*)

procedure GETEQX(var TEQX: double);
  begin
    WRITE  (' equinox (yyyy) ? ');
    READLN (TEQX); TEQX := (TEQX-2000.0)/100.0;
  end;

(*-----------------------------------------------------------------------*)

procedure GETDAT(var T: double);
  var D,M,Y  : integer;
      HOUR,JD: Double;
  begin
    WRITE (' Date (year month day hour) ?      ');
    READLN (Y,M,D,HOUR);
    JD := MJD(D,M,Y,HOUR) + 2400000.5;
    writeln; writeln (' JD',JD:13:4); writeln;
    T := (JD-2451545.0) / 36525.0;
  end;

(*-----------------------------------------------------------------------*)

procedure GETINP (var X,Y,Z,TEQX: double;var ECLIPT:Boolean);

  var I,D,M  : integer;
      L,B,R,S: Double;

  begin (* GETINP *)

    writeln;
    writeln('                COCO: coordinate conversion           ');
    writeln('                      version 93/07/01                ');
    writeln('        (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
    writeln;
    writeln (' Coordinate input: please select format required ');
    writeln;
    writeln ('  1  ecliptic   cartesian     2  ecliptic   polar');
    writeln ('  3  equatorial cartesian     4  equatorial polar');
    writeln;
    WRITE   ('  '); READLN (I); writeln;

    case I of
      1: begin
           WRITE (' Coordinates (x y z) ?  '); READLN(X,Y,Z); ECLIPT:=True;
         end;
      2: begin
           WRITE (' Coordinates (L (o '' ")  B (o '' ")  R) ?  ');
           READ(D,M,S); DDD(D,M,S,L); READLN(D,M,S,R); DDD(D,M,S,B);
           CART(R,B,L,X,Y,Z); ECLIPT:=True;
         end;
      3: begin
           WRITE (' Coordinates (x y z) ?  '); READLN(X,Y,Z);ECLIPT:=False;
         end;
      4: begin
           WRITE (' Coordinates (RA (h m s)  DEC (o '' ")  R) ?  ');
           READ(D,M,S); DDD(D,M,S,L); READLN(D,M,S,R); DDD(D,M,S,B);
           L:=L*15.0; CART(R,B,L,X,Y,Z); ECLIPT:=False;
         end;
      end; (* case *)

    GETEQX (TEQX); (* read equinox *)

  end; (* GETINP *)

(*-----------------------------------------------------------------------*)

procedure RESULT (X,Y,Z: Double; ECLIPT: Boolean);

  var L,B,R,S: Double;
      D,M    : integer;

  begin (* RESULT *)

    writeln; writeln (' (x,y,z) = (',x:13:8,',',y:13:8,',',z:13:8,')');
    writeln;

    POLAR (X,Y,Z,R,B,L);
    if ECLIPT
      then
        begin
          writeln (' ':5,'   o  ''  " ',' ':8,'   o  ''  " ');
          DMS(L,D,M,S); WRITE(' L = ',D:3,M:3,S:5:1,' ':3);
          DMS(B,D,M,S); WRITE(' B = ',D:3,M:3,S:5:1,' ':3);
        end
      else
        begin
          writeln (' ':5,'   h  m  s ',' ':10,'   o  ''  " ');
          DMS(L/15,D,M,S); WRITE(' RA = ',D:2,M:3,S:5:1,' ':3);
          DMS(B,D,M,S);    WRITE(' DEC = ',D:3,M:3,S:5:1,' ':3);
        end;

    writeln (' R = ',R:12:8); writeln; writeln;

  end; (* RESULT *)

(*-----------------------------------------------------------------------*)

begin (* COCO *)

  GETINP (X,Y,Z,TEQX,ECLIPT);
  RESULT (X,Y,Z,ECLIPT);

  repeat

    WRITE (' Command (?=Help): ');  READLN (MODE);   writeln;

    if MODE IN ['?','P','p','E','e','H','h','G','g'] then
    case MODE of

      '?'    : begin (* help *)
                 writeln;
                 writeln ('   E: equatorial <-> ecliptic                ');
                 writeln ('   P: -> precession       G: -> geocentric   ');
                 writeln ('   H: -> heliocentric     S: -> STOP         ');
                 writeln;
               end;

      'P','p': begin (* precession *)
                 WRITE (' New');
                 GETEQX(TEQXN); (* read new equinox *)
                 if ECLIPT then PMATECL(TEQX,TEQXN,A)
                           else PMATEQU(TEQX,TEQXN,A);
                 PRECART(A,X,Y,Z);
                 TEQX := TEQXN;
                 writeln;
                 writeln (' Coordinates referred to equinox T =',
                          TEQX:13:10);
               end;

      'E','e': begin (* ecliptic <-> equatorial *)
                 writeln;
                 if (ECLIPT) then
                   begin ECLEQU(TEQX,X,Y,Z); WRITE(' Equatorial'); end
                 else
                   begin EQUECL(TEQX,X,Y,Z); WRITE(' Ecliptic'); end;
                 writeln (' coordinates: ');
                 ECLIPT := NOT ECLIPT;
               end;

      'G','g', (* -> geocentric coordinates   *)
      'H','h': (* -> heliocentric coordinates *)
               begin
                 GETDAT(T); (* read date *)
                 Sun200(T,LS,BS,RS);
                 CART(RS,BS,LS,XS,YS,ZS);
                 PMATECL(T,TEQX,A);
                 PRECART(A,XS,YS,ZS);
                 if NOT ECLIPT then ECLEQU(TEQX,XS,YS,ZS);
                 if MODE IN ['G','g']
                   then (* -> geocentric *)
                     begin
                       X:=X+XS; Y:=Y+YS; Z:=Z+ZS;
                       writeln(' Geocentric coordinates: ');
                     end
                   else (* -> heliocentric *)
                     begin
                       X:=X-XS; Y:=Y-YS; Z:=Z-ZS;
                       writeln(' Heliocentric coordinates: ');
                     end;
               end;
       end;

    if NOT (MODE IN ['?','S','s']) then RESULT(X,Y,Z,ECLIPT);

  until MODE IN ['S','s'];

end. (* COCO *)


