(*-----------------------------------------------------------------------*)
(*                                 LUNA                                  *)
(*                           lunar ephemeris                             *)
(*-----------------------------------------------------------------------*)

program LUNA(Input,Output);

{$APPTYPE CONSOLE}

  uses 
   Astro.Matlib,
   Astro.Pnulib,
   Astro.Sphlib,
   Astro.Moon,
   Astro.Timlib;

  const T_OVERLAP  = 3.42E-6;        (* 3h  in Julian centuries *)
        T_DEVELOP  = 2.737850787E-4; (* 10d in Julian centuries *)

  var RA,DE,R,PAR, MODJD,HOUR : Double;
      T,DT,T_START,T_END,TA,TB: Double;
      DAY,MONTH,YEAR,NLINE    : integer;
      RA_POLY,DE_POLY,R_POLY  : TPolynom;


(*-----------------------------------------------------------------------*)
(* GETEPH: reads the period of time for the ephemeris                    *)
(*-----------------------------------------------------------------------*)
procedure GETEPH(var T1,DT,T2: double);
  var YEAR,MONTH,DAY: integer;
      HOUR          : Double;
  begin
    writeln;
    writeln('                          LUNA: lunar ephemeris            ');
    writeln('                            Version 93/07/01               ');
    writeln('             (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
    writeln;
    writeln(' Begin and end of the ephemeris: ');
    writeln;
    write  ('  first date (yyyy mm dd hh.hhh)  ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T1 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    write  ('  last date  (yyyy mm dd hh.hhh)  ');
    READLN (YEAR,MONTH,DAY,HOUR);
    T2 :=  ( MJD(DAY,MONTH,YEAR,HOUR) - 51544.5 ) / 36525.0;
    write  ('  step size (dd hh.hh)            ');
    READLN (DAY,HOUR);
    DT :=  ( DAY + HOUR/24.0 ) / 36525.0;
  end;
(*-----------------------------------------------------------------------*)
(* WRTLBRP: formatted output                                             *)
(*-----------------------------------------------------------------------*)
procedure WRTLBRP (L,B,R,P: double);
  var H,M: integer;
      S  : Double;
  begin
    DMS(L,H,M,S);  write (H:5,M:3,S:5:1);
    DMS(B,H,M,S);  write (H:5,M:3,S:5:1);  write (R:10:3);
    DMS(P,H,M,S);  if (H>0) then M:=M+60;  writeln (M:6,S:6:2);
  end;
(*-----------------------------------------------------------------------*)

begin   (* main program *)

  GETEPH(T_START,DT,T_END);  (* read desired dates *)

  writeln;
  write   ('    Date      ET       RA           Dec      Distance  ');
  writeln (' Parallax ');
  write   ('              h      h  m  s      o  ''  "  Earth radii');
  writeln ('    ''  "   ');

  T := T_START;   TB := T_START;   NLINE := 0;

  while (T<=T_END) do

    begin

      if (T>TB-T_OVERLAP) then       (* new expansion of the coordinates *)
        begin
          TA := T-T_OVERLAP;  TB := T+T_DEVELOP+T_OVERLAP;
          T_FIT_MOON (TA,TB,MAX_TP_DEG,RA_POLY,DE_POLY,R_POLY);
        end;

      (* date *)
      MODJD := T*36525.0 + 51544.5;  CALDAT (MODJD,DAY,MONTH,YEAR,HOUR);
      write  (YEAR:5,'/',MONTH:2,'/',DAY:2,HOUR:5:1);

      (* coordinates *)
      RA  := T_EVAL(RA_POLY,T)/15.0;  if RA<0.0 then RA := RA + 24.0;
      DE  := T_EVAL(DE_POLY,T);
      R   := T_EVAL(R_POLY, T);
      PAR := ASN(1.0/R);

      (* print coordinates *)
      WRTLBRP (RA,DE,R,PAR);
      NLINE := NLINE + 1;
      if (NLINE MOD 5) = 0 then writeln;

      T := T + DT;

    end;

end.

