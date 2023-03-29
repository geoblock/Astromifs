(* ----------------------------------------------------------------------- *)
(* Luna *)
(* lunar ephemeris *)
(* ----------------------------------------------------------------------- *)

program Luna(Input, Output);

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric,
  Apc.Moon,
  Apc.Time;

const
  T_OVERLAP = 3.42E-6; (* 3h  in Julian centuries *)
  T_DEVELOP = 2.737850787E-4; (* 10d in Julian centuries *)

var
  Ra, DE, R, PAR, MODJD, HOUR: Double;
  T, DT, T_START, T_END, TA, TB: Double;
  DAY, MONTH, YEAR, NLINE: integer;
  RA_POLY, DE_POLY, R_POLY: TPolynomCheb;

  (* ----------------------------------------------------------------------- *)
  (* GetEph: reads the period of time for the ephemeris *)
  (* ----------------------------------------------------------------------- *)
procedure GetEph(var T1, DT, T2: Double);
var
  YEAR, MONTH, DAY: integer;
  HOUR: Double;
begin
  writeln;
  writeln('                          Luna: lunar ephemeris            ');
  writeln('                            Version 93/07/01               ');
  writeln('             (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
  writeln;
  writeln(' Begin and end of the ephemeris: ');
  writeln;
  write('  first date (yyyy mm dd hh.hhh)  ');
  Readln(YEAR, MONTH, DAY, HOUR);
  T1 := (MJD(DAY, MONTH, YEAR, HOUR) - 51544.5) / 36525.0;
  write('  last date  (yyyy mm dd hh.hhh)  ');
  Readln(YEAR, MONTH, DAY, HOUR);
  T2 := (MJD(DAY, MONTH, YEAR, HOUR) - 51544.5) / 36525.0;
  write('  step size (dd hh.hh)            ');
  Readln(DAY, HOUR);
  DT := (DAY + HOUR / 24.0) / 36525.0;
end;

(* ----------------------------------------------------------------------- *)
(* WRTLBRP: formatted output *)
(* ----------------------------------------------------------------------- *)
procedure WRTLBRP(L, B, R, P: Double);
var
  H, M: integer;
  S: Double;
begin
  DMS(L, H, M, S);
  write(H:5, M:3, S:5:1);
  DMS(B, H, M, S);
  write(H:5, M:3, S:5:1);
  write(R:10:3);
  DMS(P, H, M, S);
  if (H > 0) then
    M := M + 60;
  writeln(M:6, S:6:2);
end;
(* ----------------------------------------------------------------------- *)

begin (* main program *)

  GetEph(T_START, DT, T_END); (* read desired dates *)

  writeln;
  write('    Date      ET       Ra           Dec      Distance  ');
  writeln(' Parallax ');
  write('              h      h  m  s      o  ''  "  Earth radii');
  writeln('    ''  "   ');

  T := T_START;
  TB := T_START;
  NLINE := 0;

  while (T <= T_END) do

  begin

    if (T > TB - T_OVERLAP) then (* new expansion of the coordinates *)
    begin
      TA := T - T_OVERLAP;
      TB := T + T_DEVELOP + T_OVERLAP;
      T_Fit_Moon(TA, TB, MAX_TP_DEG, RA_POLY, DE_POLY, R_POLY);
    end;

    (* date *)
    MODJD := T * 36525.0 + 51544.5;
    CalDat(MODJD, DAY, MONTH, YEAR, HOUR);
    write(YEAR:5, '/', MONTH:2, '/', DAY:2, HOUR:5:1);

    (* coordinates *)
    Ra := T_Eval(RA_POLY, T) / 15.0;
    if Ra < 0.0 then
      Ra := Ra + 24.0;
    DE := T_Eval(DE_POLY, T);
    R := T_Eval(R_POLY, T);
    PAR := ASN(1.0 / R);

    (* print coordinates *)
    WRTLBRP(Ra, DE, R, PAR);
    NLINE := NLINE + 1;
    if (NLINE MOD 5) = 0 then
      writeln;

    T := T + DT;

  end;

end.
