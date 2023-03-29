program Coco(Input, Output);

//-----------------------------------------------------------------------
(* Coordinate conversion *)
//-----------------------------------------------------------------------

{$APPTYPE CONSOLE}

uses
  Apc.Mathem,
  Apc.PrecNut,
  Apc.Spheric,
  Apc.Sun,
  Apc.Time;

var
  X, Y, Z, XS, YS, ZS: Double;
  T, TEQX, TEQXN: Double;
  LS, BS, RS: Double;
  A: Double33;
  IsEclipt: Boolean;
  CharMode: Char;

//-----------------------------------------------------------------------

procedure GetEqx(var TEQX: Double);
begin
  write(' equinox (yyyy) ? ');
  Readln(TEQX);
  TEQX := (TEQX - 2000.0) / 100.0;
end;

//-----------------------------------------------------------------------

procedure GetDat(var T: Double);
var
  D, M, Y: integer;
  HOUR, JD: Double;
begin
  write(' Date (year month day hour) ?      ');
  Readln(Y, M, D, HOUR);
  JD := MJD(D, M, Y, HOUR) + 2400000.5;
  writeln;
  writeln(' JD', JD:13:4);
  writeln;
  T := (JD - 2451545.0) / 36525.0;
end;

//-----------------------------------------------------------------------

procedure GetInp(var X, Y, Z, TEQX: Double; var IsEclipt: Boolean);
var
  I, D, M: integer;
  L, B, R, S: Double;
begin

  writeln;
  writeln('                Coco: coordinate conversion           ');
  writeln('                      version 93/07/01                ');
  writeln('        (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
  writeln;
  writeln(' Coordinate input: please select format required ');
  writeln;
  writeln('  1  ecliptic   cartesian     2  ecliptic   polar');
  writeln('  3  equatorial cartesian     4  equatorial polar');
  writeln;
  write('  ');
  Readln(I);
  writeln;

  case I of
    1:
      begin
        write(' Coordinates (x y z) ?  ');
        Readln(X, Y, Z);
        IsEclipt := True;
      end;
    2:
      begin
        write(' Coordinates (L (o '' ")  B (o '' ")  R) ?  ');
        read(D, M, S);
        Ddd(D, M, S, L);
        Readln(D, M, S, R);
        Ddd(D, M, S, B);
        Cart(R, B, L, X, Y, Z);
        IsEclipt := True;
      end;
    3:
      begin
        write(' Coordinates (x y z) ?  ');
        Readln(X, Y, Z);
        IsEclipt := False;
      end;
    4:
      begin
        write(' Coordinates (Ra (h m s)  Dec (o '' ")  R) ?  ');
        read(D, M, S);
        Ddd(D, M, S, L);
        Readln(D, M, S, R);
        Ddd(D, M, S, B);
        L := L * 15.0;
        Cart(R, B, L, X, Y, Z);
        IsEclipt := False;
      end;
  end; // case
  GetEqx(TEQX); // read equinox
end; // GetInp

//-----------------------------------------------------------------------

procedure Results(X, Y, Z: Double; IsEclipt: Boolean);

var
  L, B, R, S: Double;
  D, M: integer;

begin (* Results *)

  writeln;
  writeln(' (x,y,z) = (', X:13:8, ',', Y:13:8, ',', Z:13:8, ')');
  writeln;

  Polar(X, Y, Z, R, B, L);
  if IsEclipt then
  begin
    writeln(' ':5, '   o  ''  " ', ' ':8, '   o  ''  " ');
    DMS(L, D, M, S);
    write(' L = ', D:3, M:3, S:5:1, ' ':3);
    DMS(B, D, M, S);
    write(' B = ', D:3, M:3, S:5:1, ' ':3);
  end
  else
  begin
    writeln(' ':5, '   h  m  s ', ' ':10, '   o  ''  " ');
    DMS(L / 15, D, M, S);
    write(' Ra = ', D:2, M:3, S:5:1, ' ':3);
    DMS(B, D, M, S);
    write(' Dec = ', D:3, M:3, S:5:1, ' ':3);
  end;

  writeln(' R = ', R:12:8);
  writeln;
  writeln;
end; // Results

//-----------------------------------------------------------------------

begin (* Coco *)
  GetInp(X, Y, Z, TEQX, IsEclipt);
  Results(X, Y, Z, IsEclipt);
  repeat
    write(' Command (?=Help): ');
    Readln(CharMode);
    writeln;

    if CharMode in ['?', 'P', 'p', 'E', 'e', 'H', 'h', 'G', 'g'] then
      case CharMode of
        '?':
          begin (* help *)
            writeln;
            writeln('   E: equatorial <-> ecliptic                ');
            writeln('   P: -> precession       G: -> geocentric   ');
            writeln('   H: -> heliocentric     S: -> STOP         ');
            writeln;
          end;
        'P', 'p':
          begin (* precession *)
            write(' New');
            GetEqx(TEQXN); (* read new equinox *)
            if IsEclipt then
              PrecMatEcl(TEQX, TEQXN, A)
            else
              PrecMatEqu(TEQX, TEQXN, A);
            PrecArt(A, X, Y, Z);
            TEQX := TEQXN;
            writeln;
            writeln(' Coordinates referred to equinox T =', TEQX:13:10);
          end;
        'E', 'e':
          begin (* ecliptic <-> equatorial *)
            writeln;
            if (IsEclipt) then
            begin
              Ecl2Equ(TEQX, X, Y, Z);
              write(' Equatorial');
            end
            else
            begin
              Equ2Ecl(TEQX, X, Y, Z);
              write(' Ecliptic');
            end;
            writeln(' coordinates: ');
            IsEclipt := not IsEclipt;
          end;
        'G', 'g', (* -> geocentric coordinates *)
        'H', 'h': (* -> heliocentric coordinates *)
          begin
            GetDat(T); (* read date *)
            SunPos(T, LS, BS, RS);
            Cart(RS, BS, LS, XS, YS, ZS);
            PrecMatEcl(T, TEQX, A);
            PrecArt(A, XS, YS, ZS);
            if not IsEclipt then
              Ecl2Equ(TEQX, XS, YS, ZS);
            if CharMode in ['G', 'g'] then (* -> geocentric *)
            begin
              X := X + XS;
              Y := Y + YS;
              Z := Z + ZS;
              writeln(' Geocentric coordinates: ');
            end
            else (* -> heliocentric *)
            begin
              X := X - XS;
              Y := Y - YS;
              Z := Z - ZS;
              writeln(' Heliocentric coordinates: ');
            end;
          end;
      end;
    if not(CharMode in ['?', 'S', 's']) then
      Results(X, Y, Z, IsEclipt);
  until CharMode in ['S', 's'];
end.
