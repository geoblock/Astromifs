(*-----------------------------------------------------------------------*)
(*                                ECLIPSE                                *)
(*             central line and duration of a solar eclipse              *)
(*-----------------------------------------------------------------------*)
program ECLIPSE(INPUT,OUTPUT);

  uses 
    Astro.Matlib, Astro.Pnulib, Astro.Sphlib, Astro.Sunlib, Astro.Moonlib, Astro.Timlib;

  const MAX_TP_DEG = 8;               (* max. degree of Chebyshev polyn. *)
        H          = 1.14E-6;         (* 1h in jul.cent. (1/(24*36525))  *)

  type  PHASE_TYPE = ( NO_ECLIPSE, PARTIAL,
                       NON_CEN_ANN, NON_CEN_TOT, ANNULAR, TOTAL );

  var   T_begin,T_END,T,DT,MJDUT             : Double;
        ETDIFUT                              : Double;
        LAMBDA,PHI,T_UMBRA                   : Double;
        RAM_POLY,DEM_POLY,RM_POLY            : TPolynom;
        RAS_POLY,DES_POLY,RS_POLY            : TPolynom;
        PHASE                                : PHASE_TYPE;


(*-----------------------------------------------------------------------*)
(* GET_INPUT: determines search interval, step size and ET-UT            *)
(* (T_begin,T_END,DT in Julian centuries since J2000 UT; ET-UT in sec)   *)
(*-----------------------------------------------------------------------*)

procedure GET_INPUT ( var T_begin, T_END, DT, ETDIFUT: Double);

  var D,M,Y     : integer;
      UT,T,DTAB : Double;
      VALID     : Boolean;

  begin

    writeln;
    writeln ('                     ECLIPSE: solar eclipses            ');
    writeln ('                            93/07/01                    ');
    writeln ('          (c) 1993 Thomas Pfleger, Oliver Montenbruck   ');
    writeln;
    WRITE   (' Date of New Moon (yyyy mm dd UT): '); READLN(Y,M,D,UT);
    WRITE   (' Output step size (min)          : '); READLN(DTAB);
    DT := (DTAB/1440.0)/36525.0;
    UT := TRUNC ( UT*60.0/DTAB + 0.5 ) * DTAB/60.0 ; (* round to 1 min *)
    T := (MJD(D,M,Y,UT)-51544.5)/36525.0;
    T_begin := T-0.25/36525.0;   T_END := T+0.25/36525.0;
    ETMINUT ( T, ETDIFUT, VALID);
    WRITE   (' Difference ET-UT (sec)          : ');
    if (VALID) then  WRITE ('  (proposal:',Trunc(ETDIFUT+0.5):4,' sec) ');
    READLN(ETDIFUT);
    writeln;
    writeln ('    Date       UT         Phi     Lambda   Durat  Phase ');
    writeln;
    writeln ('               h  m        o  ''      o  ''   min       ');

  end;


(*-----------------------------------------------------------------------*)
(* WRTOUT: formated output                                               *)
(*-----------------------------------------------------------------------*)

procedure WRTOUT ( MJDUT,LAMBDA,PHI,T_UMBRA: Double; PHASE: PHASE_TYPE );

  var DAY,MONTH,YEAR,H,M: integer;
      HOUR,S            : Double;

  begin

    CALDAT ( MJDUT,DAY,MONTH,YEAR,HOUR );      (* date  *)
    WRITE  ( YEAR:5,'/',MONTH:2,'/',DAY:2);
    DMS(HOUR+0.5/60.0,H,M,S); WRITE (H:4,M:3); (* time rounded to 1 min *)
    if ( ORD(PHASE)<ORD(ANNULAR) )
      then
        WRITE('       -- --     -- --    ---')
      else
        begin
          DMS(PHI,H,M,S);    WRITE (H:9,M:3);
          DMS(LAMBDA,H,M,S); WRITE (H:7,M:3); WRITE (T_UMBRA:7:1);
        end;
    case PHASE of
      NO_ECLIPSE : WRITE('   ----                  ');
      PARTIAL    : WRITE('   partial               ');
      NON_CEN_ANN: WRITE('   annular (non-central) ');
      NON_CEN_TOT: WRITE('   total (non-central)   ');
      ANNULAR    : WRITE('   annular               ');
      TOTAL      : WRITE('   total                 ');
    end;
    writeln;
  end;

(*-----------------------------------------------------------------------*)
(* INTSECT: calculates the intersection of the shadow axis with the      *)
(*          surface of the Earth                                         *)
(*                                                                       *)
(*   RAM,DEM,RM,  equatorial coordinates of Moon and Sun (right asc. RA  *)
(*   RAS,DES,RS:  and declination in deg; distance in Earth radii)       *)
(*   X,Y,Z:       equatorial coord. of the shadow point (in Earth radii) *)
(*   EX,EY,EZ:    unit vector of the shadow axis                         *)
(*   D_UMBRA:     umbra diameter in Earth radii                          *)
(*   PHASE:       phase of the eclipse                                   *)
(*-----------------------------------------------------------------------*)

procedure INTSECT ( RAM,DEM,RM, RAS,DES,RS: Double;
                    var X,Y,Z, EX,EY,EZ, D_UMBRA: Double;
                    var PHASE: PHASE_TYPE );

  const FAC = 0.996633;     (* ratio polar/equatorial Earth radius       *)
        D_M =   0.5450;     (* lunar diameter in units of 1 Earth radius *)
        D_S = 218.25;       (* solar diameter in units of 1 Earth radius *)

  var XM,YM,ZM, XS,YS,ZS, XMS,YMS,ZMS, RMS: Double;
      DELTA, R0, S0, S, D_PENUMBRA        : Double;

  begin

    CART(RM,DEM,RAM,XM,YM,ZM); ZM:=ZM/FAC; (* solar and lunar coordinat. *)
    CART(RS,DES,RAS,XS,YS,ZS); ZS:=ZS/FAC; (* scale z-coordinate         *)

    XMS:=XM-XS;  YMS:=YM-YS;  ZMS:=ZM-ZS;   (* vector Sun -> Moon,       *)
    RMS := Sqrt(XMS*XMS+YMS*YMS+ZMS*ZMS);   (* distance Sun -> Moon      *)
    EX:=XMS/RMS;  EY:=YMS/RMS;  EZ:=ZMS/RMS;(* unit vector Sun -> Moon   *)

    S0 := -( XM*EX + YM*EY + ZM*EZ ); (* dist. Moon -> fundamental plane *)
    DELTA := S0*S0+1.0-XM*XM-YM*YM-ZM*ZM;
    R0 := Sqrt(1.0-DELTA);            (* dist.center Earth - shadow axis *)

    D_UMBRA    := (D_S-D_M)*S0/RMS-D_M;     (* diameter of pen-/umbra    *)
    D_PENUMBRA := (D_S+D_M)*S0/RMS+D_M;     (* on the fundamental plane  *)


    (* determine phase and shadow coordinates if required                *)
    if ( R0 < 1.0 )
      then                       (* shadow axis intersects the Earth     *)
        begin                    (* -> total or annular eclipse          *)
          S := S0-Sqrt(DELTA);
          D_UMBRA := (D_S-D_M)*(S/RMS)-D_M;   (* umbra diameter on Earth *)
          X:=XM+S*EX;  Y:=YM+S*EY; Z:=ZM+S*EZ;
          Z:=Z*FAC;                           (* rescale z-coordinate    *)
          if D_UMBRA>0 then PHASE:=ANNULAR else PHASE:=TOTAL;
        end
      else
        if ( R0 < 1.0+0.5*abs(D_UMBRA) )
          then                                (* non-central eclipse     *)
            if D_UMBRA>0 then PHASE:=NON_CEN_ANN else PHASE:=NON_CEN_TOT
          else
            if ( R0 < 1.0+0.5*D_PENUMBRA)
              then PHASE := PARTIAL           (* partial eclipse         *)
              else PHASE := NO_ECLIPSE;       (* no eclipse              *)


  end;


(*-----------------------------------------------------------------------*)
(* CENTRAL: central line, phase and duration of the eclipse              *)
(*                                                                       *)
(*   T_UT:        time in Julian centuries since J2000 UT                *)
(*   ETDIFUT:     difference Ephemeris Time - Universal Time (in sec)    *)
(*   RAM_POLY,DEM_POLY,RM_POLY, RAS_POLY,DES_POLY,RS_POLY:               *)
(*                Chebyshev coefficients for solar and lunar coordinates *)
(*   LAMBDA, PHI: geographic long. and latit. of the shadow center (deg) *)
(*   T_UMBRA:     duration of the total or annular phase (min)           *)
(*   PHASE:       phase of the eclipse                                   *)
(*-----------------------------------------------------------------------*)

procedure CENTRAL ( T_UT, ETDIFUT            : Double;
                    RAM_POLY,DEM_POLY,RM_POLY: TPolynom;
                    RAS_POLY,DES_POLY,RS_POLY: TPolynom;
                    var LAMBDA,PHI,T_UMBRA   : Double;
                    var PHASE                : PHASE_TYPE );

  const AU    = 23454.78;     (* 1AU in earth radii (149597870/6378.14)  *)
        DT    = 0.1;          (* small time interval; dt = 0.1 min       *)
        MPC   = 52596000.0;   (* minutes per Julian century (1440*36525) *)
        OMEGA = 4.3755E-3;    (* angular velocity of the earth (rad/min) *)

  var RAM,DEM,RM, RAS,DES,RS, RA,DEC,R, DX,DY,DZ,D, MJDUT  : Double;
      T,X,Y,Z,EX,EY,EZ,D_UMBRA, XX,YY,ZZ,EXX,EYY,EZZ,DU, W : Double;
      PH                                                   : PHASE_TYPE;

  (* calculate lunar and solar coordinates from Chebyshev coefficients   *)
  procedure POSITION ( T: Double; var RAM,DEM,RM, RAS,DES,RS : Double);
    begin
      RAM:=T_EVAL(RAM_POLY,T); RAS:=T_EVAL(RAS_POLY,T);
      DEM:=T_EVAL(DEM_POLY,T); DES:=T_EVAL(DES_POLY,T);
      RM :=T_EVAL(RM_POLY,T);  RS :=T_EVAL(RS_POLY,T);  RS:=RS*AU;
    end;

  begin

    (* julian centuries since J2000 ET *)
    T := T_UT + ETDIFUT/(86400.0*36525.0);

    (* phase of eclipse and coordinates of the shadow at time T *)
    POSITION ( T, RAM,DEM,RM, RAS,DES,RS );
    INTSECT ( RAM,DEM,RM, RAS,DES,RS, X,Y,Z,EX,EY,EZ, D_UMBRA,PHASE );

    (* for central phase only:  geogr. coord. and duration of totality *)
    if ( ORD(PHASE) < ORD(ANNULAR) )
      then begin LAMBDA:=0.0; PHI:=0.0; T_UMBRA:=0.0; end
      else
        begin
          (* geographic coordinates: *)
          MJDUT := 36525.0*T_UT + 51544.5;
          POLAR ( X,Y,Z, R,DEC,RA );
          PHI    := DEC + 0.1924*SN(2.0*DEC);
          LAMBDA := 15.0*LMST(MJDUT,0.0)-RA;
          if LAMBDA>+180.0 then LAMBDA:=LAMBDA-360.0;
          if LAMBDA<-180.0 then LAMBDA:=LAMBDA+360.0;
          (* duration of totality for this place           *)
          (* (a) shadow coordinates at time T+DT (or T-DT) *)
          POSITION ( T+DT/MPC, RAM,DEM,RM, RAS,DES,RS );  W:=+DT*OMEGA;
          INTSECT ( RAM,DEM,RM, RAS,DES,RS, XX,YY,ZZ,EXX,EYY,EZZ, DU, PH );
          if (ORD(PH)<ORD(ANNULAR)) then
            begin
              POSITION ( T-DT/MPC, RAM,DEM,RM,RAS,DES,RS);  W:=-DT*OMEGA;
              INTSECT (RAM,DEM,RM,RAS,DES,RS, XX,YY,ZZ,EXX,EYY,EZZ, DU,PH);
            end;
          (* (b) displacement DX,DY,DZ of the shadow on Earth    *)
          (*     and fraction D perpendicular to the shadow axis *)
          DX := XX-X+W*Y;   DY := YY-Y-W*X;   DZ := ZZ-Z;
          D  := Sqrt( DX*DX+DY*DY+DZ*DZ -
                      (DX*EX+DY*EY+DZ*EZ)*(DX*EX+DY*EY+DZ*EZ) );
          T_UMBRA := DT * abs(D_UMBRA) / D;
        end;

  end;

(*-----------------------------------------------------------------------*)

begin (* main program *)


  (* read search interval *)

  GET_INPUT ( T_begin,T_END, DT, ETDIFUT );


  (* Chebyshev approximations *)

  T_FIT_MOON (T_begin-H,T_END+H,8,RAM_POLY,DEM_POLY,RM_POLY);
  T_FIT_SUN  (T_begin-H,T_END+H,3,RAS_POLY,DES_POLY,RS_POLY);


  (* calculate phase and central line of the eclipse *)

  T := T_begin;

  repeat

    CENTRAL ( T,ETDIFUT,
              RAM_POLY,DEM_POLY,RM_POLY, RAS_POLY,DES_POLY,RS_POLY,
              LAMBDA, PHI, T_UMBRA, PHASE  );

    if PHASE<>NO_ECLIPSE then
      begin
        MJDUT := 36525.0*T + 51544.5;
        WRTOUT ( MJDUT, LAMBDA, PHI, T_UMBRA, PHASE );
      end;

    T := T+DT;

  until (T > T_END);

end.
