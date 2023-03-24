(*-----------------------------------------------------------------------*)
(*                              PLANPOS                                  *)
(*           heliocentric and geocentric planetary positions             *)
(*-----------------------------------------------------------------------*)

program PLANPOS(INPUT,OUTPUT);

  uses 
    Astro.Matlib, Astro.Pnulib, Astro.Sphlib, Astro.Sunlib, 
    Astro.p15lib, Astro.p69lib, Astro.Timlib;

  const J2000 =  0.0;
        B1950 = -0.500002108;

  var DAY, MONTH, YEAR, IPLAN, IMODE, K    : integer;
      HOUR, MODJD, T, TEQX                 : Double;
      X,Y,Z, XP,YP,ZP, XS,YS,ZS            : Double;
      L,B,R, LS,BS,RS, RA,DEC,DELTA,DELTA0 : Double;
      A                                    : Double33;
      MODE                                 : CHAR;

(*-----------------------------------------------------------------------*)
(* PRINTOUT: print coordinates of one planet                             *)
(*-----------------------------------------------------------------------*)
procedure PRINTOUT(IPLAN:Integer;L,B,R,RA,DEC,DELTA: double);
  var H,M: integer;
      S  : Double;
  begin
    DMS(L,H,M,S);   WRITE (H:3,M:3,S:5:1);
    DMS(B,H,M,S);   WRITE (H:4,M:3,S:5:1);
    if IPLAN<4 then WRITE (R:11:6)
               else WRITE (R:10:5,' ');
    DMS(RA,H,M,S);  WRITE (H:4,M:3,S:6:2);
    DMS(DEC,H,M,S); WRITE (H:4,M:3,S:5:1);
    if IPLAN<4 then WRITE (DELTA:11:6)
               else WRITE (DELTA:10:5,' ');
  end;
(*-----------------------------------------------------------------------*)

begin (* PLANPOS *)

  writeln;
  writeln('   PLANPOS: geocentric and heliocentric planetary positions ');
  writeln('                       version 93/07/01                     ');
  writeln('         (c) 1993 Thomas Pfleger, Oliver Montenbruck        ');
  writeln;

  repeat

    writeln;
    writeln (' (J) J2000 astrometric       (B) B1950 astrometric  ');
    writeln (' (A) apparent coordinates    (E) end                ');
    writeln;
    WRITE   (' enter option: '); READLN (MODE);
    writeln;

    if MODE IN ['A','a','J','j','B','b']  then

      begin

        (* read date *)

        WRITE (' date (year month day hour) ?     ');
        READLN (YEAR,MONTH,DAY,HOUR);  writeln; writeln; writeln; writeln;
        MODJD := MJD (DAY,MONTH,YEAR,HOUR);  T:=(MODJD-51544.5)/36525.0;
        WRITE (' date:  ', YEAR:4,'/',MONTH:2,'/',DAY:2,' ',HOUR:5:1,'(ET)');
        WRITE ('JD:':6,(MODJD+2400000.5):12:3,'equinox ':18);
        case MODE of
          'A','a': writeln ('of date');
          'J','j': writeln ('J2000');
          'B','b': writeln ('B1950');
         end;
        writeln;

        (* header *)

        WRITE   (' ':10,'l':6,'b':12,'r':11);
        writeln (' ':7,'RA':5,'Dec':13,'delta':13);
        WRITE   (' ':9,'   o  ''  "',' ':3,'  o  ''  "',' ':6,'AU',' ':4);
        writeln (' ':2,'  h  m  s',' ':4,'  o  ''  "',' ':6,'AU');

        (* ecliptic coordinates of the sun, Equinox T  *)

        Sun200 (T,LS,BS,RS);

        (* planetary coordinates; include Pluto between 1890 and 2100  *)

        if ( (-1.1<T) and (T<+1.0) ) then K:=9 else K:=8;

        for IPLAN:=0 to K do

          begin

            (* heliocentric ecliptic coordinates of the planet        *)

            case IPLAN of
              1: Mercury200(T,L,B,R); 2: Venus200(T,L,B,R);
              4: Mars200(T,L,B,R); 5: Jupiter200(T,L,B,R); 6: Saturn200(T,L,B,R);
              7: Uranus200(T,L,B,R); 8: Neptune200(T,L,B,R); 9: Pluto200(T,L,B,R);
              0: begin L:=0.0; B:=0.0; R:=0.0; end;
              3: begin L:=LS+180.0; B:=-BS; R:=RS; end;
             end;

            (* geocentric ecliptic coordinates (light-time corrected)  *)

            if MODE IN ['A','a'] then IMODE:=2 else IMODE:=1;
            GEOCEN (T, L,B,R, LS,BS,RS, IPLAN,IMODE,
                    XP,YP,ZP, XS,YS,ZS, X,Y,Z,DELTA0);

            (* precession, equatorial coordinates, nutation            *)

            case MODE of
              'J','j': TEQX:=J2000; 'B','b': TEQX:=B1950;
             end;

            if MODE IN ['A','a']
              then
                begin  ECLEQU(T,X,Y,Z); NUTEQU(T,X,Y,Z); end
              else
                begin
                  PMATECL(T,TEQX,A); PRECART (A,XP,YP,ZP);
                  PRECART (A,X,Y,Z); ECLEQU (TEQX,X,Y,Z);
                end;

            (* spherical coordinates *)

            POLAR (XP,YP,ZP,R,B,L);
            POLAR (X,Y,Z,DELTA,DEC,RA); RA:=RA/15.0;

            (* output *)

            case IPLAN of
              0: WRITE(' Sun     ');  1: WRITE(' Mercury ');
              2: WRITE(' Venus   ');  3: WRITE(' Earth   ');
              4: WRITE(' Mars    ');  5: WRITE(' Jupiter ');
              6: WRITE(' Saturn  ');  7: WRITE(' Uranus  ');
              8: WRITE(' Neptune ');  9: WRITE(' Pluto   ');
             end;

            PRINTOUT(IPLAN,L,B,R,RA,DEC,DELTA0); writeln;

          end;

          writeln;
          writeln (' l,b,r:   heliocentric ecliptic (geometric) ');
          WRITE   (' RA,Dec:  geocentric equatorial ');
          if MODE IN ['A','a'] then writeln('(apparent)')
                               else writeln('(astrometric)');
          writeln (' delta:   geocentric distance   (geometric)');
          writeln;

      end;

  until MODE IN ['E','e']

end. (* PLANPOS *)

