program Sunset(Input, Output);

(* ----------------------------------------------------------------------- *)
(* Solar and lunar rising and setting times *)
(* ----------------------------------------------------------------------- *)

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.Sun,
  Apc.Moon,
  Apc.Time;

var
  ABOVE, RISE, SETT: Boolean;
  DAY, MONTH, YEAR, I, IOBJ, NZ: Integer;
  LAMBDA, ZONE, PHI, SPHI, CPHI: Double;
  TSTART, DATE, HOUR, HH, UTRISE, UTSET: Double;
  Y_MINUS, Y_0, Y_PLUS, ZERO1, ZERO2, XE, YE: Double;
  SINH0: array [1 .. 3] of Double;

  (* ----------------------------------------------------------------------- *)
  (* SIN_ALT: sin(altitude) *)
  (* IOBJ:  1=moon, 2=sun *)
  (* ----------------------------------------------------------------------- *)
function SIN_ALT(IOBJ: integer; MJD0, HOUR, LAMBDA, CPHI, SPHI: Double): Double;
var
  MJD, T, Ra, Dec, TAU: Double;
begin
  MJD := MJD0 + HOUR / 24.0;
  T := (MJD - 51544.5) / 36525.0;
  if (IOBJ = 1) then
    MiniMoon(T, Ra, Dec)
  else
    MiniSun(T, Ra, Dec);
  TAU := 15.0 * (LMST(MJD, LAMBDA) - Ra);
  SIN_ALT := SPHI * SN(Dec) + CPHI * CS(Dec) * CS(TAU);
end;

(* ----------------------------------------------------------------------- *)
(* WHM: write time in hours and minutes *)
(* ----------------------------------------------------------------------- *)
procedure WHM(UT: Double);
var
  H, M: integer;
begin
  UT := Trunc(UT * 60.0 + 0.5) / 60.0; (* round to 1 min *)
  H := Trunc(UT);
  M := Trunc(60.0 * (UT - H) + 0.5);
  write(H:5, ':', M:2, ' ');
end;

(* ----------------------------------------------------------------------- *)
(* GetInp: read input data *)
(* ----------------------------------------------------------------------- *)
procedure GetInp(var DATE, LAMBDA, PHI, ZONE: Double);
var
  D, M, Y: integer;
begin
  writeln;
  writeln('      Sunset: solar and lunar rising and setting times ');
  writeln('                      version 93/07/01                 ');
  writeln('        (c) 1993 Thomas Pfleger, Oliver Montenbruck    ');
  writeln;
  writeln;
  write(' First date (yyyy mm dd)               ... ');
  readln(Y, M, D);
  writeln;
  write(' Observing site:  longitude (l<0 east) ... ');
  readln(LAMBDA);
  write('                  latitude             ... ');
  readln(PHI);
  write('                  local time - UT (h)  ... ');
  readln(ZONE);
  writeln;
  writeln;
  write('    Date');
  writeln('             Moon              Sun             Twilight ');
  writeln;
  write('        ');
  writeln('           rise/set          rise/set        end/beginning');
  writeln;
  ZONE := ZONE / 24.0;
  DATE := MJD(D, M, Y, 0) - ZONE;
end;

(* ----------------------------------------------------------------------- *)

begin (* Sunset *)
  SINH0[1] := SN(+8.0 / 60.0); (* moonrise          at h= +8' *)
  SINH0[2] := SN(-50.0 / 60.0); (* sunrise           at h=-50' *)
  SINH0[3] := SN(-12.0); (* nautical twilight at h=-12 degrees *)

  GetInp(TSTART, LAMBDA, PHI, ZONE);

  SPHI := SN(PHI);
  CPHI := CS(PHI);
  for I := 0 to 9 do (* loop over 10 subsequent days *)
  begin
    DATE := TSTART + I;
    CalDat(DATE + ZONE, DAY, MONTH, YEAR, HH);
    write(YEAR:5, '/', MONTH:2, '/', DAY:2, '  '); (* print current date *)
    for IOBJ := 1 to 3 do
    begin
      HOUR := 1.0;
      Y_MINUS := SIN_ALT(IOBJ, DATE, HOUR - 1.0, LAMBDA, CPHI, SPHI) - SINH0[IOBJ];

      ABOVE := (Y_MINUS > 0.0);
      RISE := False;
      SETT := False;

      (* loop over search intervals from [0h-2h] to [22h-24h] *)
      repeat
        Y_0 := SIN_ALT(IOBJ, DATE, HOUR, LAMBDA, CPHI, SPHI) - SINH0[IOBJ];
        Y_PLUS := SIN_ALT(IOBJ, DATE, HOUR + 1.0, LAMBDA, CPHI, SPHI) - SINH0[IOBJ];
        (* find parabola through three values Y_MINUS,Y_0,Y_PLUS *)
        Quad(Y_MINUS, Y_0, Y_PLUS, XE, YE, ZERO1, ZERO2, NZ);
        case (NZ) of
          0:
            ;
          1:
            if (Y_MINUS < 0.0) then
            begin
              UTRISE := HOUR + ZERO1;
              RISE := True;
            end
            else
            begin
              UTSET := HOUR + ZERO1;
              SETT := True;
            end;
          2:
            begin
              if (YE < 0.0) then
              begin
                UTRISE := HOUR + ZERO2;
                UTSET := HOUR + ZERO1;
              end
              else
              begin
                UTRISE := HOUR + ZERO1;
                UTSET := HOUR + ZERO2;
              end;
              RISE := True;
              SETT := True;
            end;
        end;

        Y_MINUS := Y_PLUS; (* prepare for next interval *)
        HOUR := HOUR + 2.0;

      until ((HOUR = 25.0) or (RISE and SETT));
      if (RISE or SETT) (* output *)
      then
      begin
        if RISE then
          WHM(UTRISE)
        else
          write('----- ':9);
        if SETT then
          WHM(UTSET)
        else
          write('----- ':9);
      end
      else
      begin
        if ABOVE then
          case IOBJ of
            1, 2:
              write('   always visible ');
            3:
              write('    always bright ');
          end
        else
          case IOBJ of
            1, 2:
              write('  always invisible');
            3:
              write('     always dark  ');
          end;
      end;
    end;
    writeln;
  end; (* end of loop over 10 days *)
  writeln;
  write(' all times in local time ( = UT ');
  if ZONE >= 0 then
    write('+');
  writeln(ZONE * 24.0:5:1, 'h )');

end.
