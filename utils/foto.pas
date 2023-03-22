(*-----------------------------------------------------------------------*)
(*                                 FOTO                                  *)
(*              astrometric analysis of photographic plates              *)
(*-----------------------------------------------------------------------*)
program FOTO (INPUT,OUTPUT,FOTINP);

  USES  
    Astro.Matlib, Astro.Sphlib;

  const MAXDIM      = 30;      (* maximum number of objects on the photo *)
        NAME_LENGTH = 12;
        ARC         = 206264.8;                 (* arcseconds per radian *)

  type  NAME_TYPE  =  array[1..NAME_LENGTH] of CHAR;
        REAL_ARRAY =  array[1..MAXDIM] of Double;
        NAME_ARRAY =  array[1..MAXDIM] of NAME_TYPE;

  var I,J,K, NREF,NOBJ, DEG,MIN  : integer;
      RA0,DEC0, A,B,C,D,E,F, SEC : Double;
      RA_OBS,DEC_OBS, D_RA,D_DEC : Double;
      DET, FOC_LEN, SCALE        : Double;
      RA,DEC, X,Y, XX,YY, DELTA  : REAL_ARRAY;
      S                          : LSQVEC;
      AA                         : LSQMAT;
      NAME                       : NAME_ARRAY;
      FOTINP                     : TEXT;


(*-----------------------------------------------------------------------*)
(* GETINP: read input data from file FOTINP                              *)
(*-----------------------------------------------------------------------*)
procedure GETINP ( var RA0,DEC0: Double;  var NOBJ: integer;
                   var NAME: NAME_ARRAY;  var RA,DEC,X,Y: REAL_ARRAY );

  var I,K,H,M: integer;
      S      : Double;
      C      : CHAR;

  begin

    writeln;
    writeln ('     FOTO: astrometric analysis of photographic plates  ');
    writeln ('                      version 93/07/01                  ');
    writeln ('        (c) 1993 Thomas Pfleger, Oliver Montenbruck     ');
    writeln;

    (* open file for reading *
    (* RESET(FOTINP); *)                              (* Standard Pascal *)
    ASSIGN(FOTINP,'FOTINP.DAT'); RESET(FOTINP);       (* TURBO Pascal    *)

    writeln (' Input data file: FOTINP'); writeln;

    (*read coordinates of the plate center *)
    for K:=1 to NAME_LENGTH do READ(FOTINP,C);
    READ  (FOTINP,H,M,S); DDD(H,M,S,RA0); RA0:=15.0*RA0;
    READLN(FOTINP,H,M,S); DDD(H,M,S,DEC0);

    (* read name, plate coordinates (and equatorial coordinates) *)
    I := 0;
    repeat
      I := I+1;
      for K:=1 to NAME_LENGTH do READ(FOTINP,NAME[I][K]);   (* name *)
      if NAME[I][1]='*'
        then (* reference star *)
          begin
            READ  (FOTINP,X[I],Y[I]);
            READ  (FOTINP,H,M,S); DDD(H,M,S,RA[I]); RA[I]:=15.0*RA[I];
            READLN(FOTINP,H,M,S); DDD(H,M,S,DEC[I]);
          end
        else (* unknown object *)
          begin
            READLN (FOTINP,X[I],Y[I]);  RA[I]:=0.0; DEC[I]:=0.0;
          end;
    until EOF(FOTINP);

    NOBJ := I;

  end;

(*-----------------------------------------------------------------------*)

begin  (* FOTO *)

  (* read input from file *)

  GETINP ( RA0,DEC0, NOBJ, NAME, RA,DEC,X,Y );

  (* calculate standard coordinates of reference stars; *)
  (* fill elements of matrix AA (column AA[*,5] serves  *)
  (* as intermediate storage)                           *)

  J:=0;
  for I:=1 to NOBJ do
    if NAME[I][1]='*' then
      begin
        J := J+1;
        EQUSTD ( RA0,DEC0, RA[I],DEC[I], XX[I],YY[I] );
        AA[J,1]:=  X[I]; AA[J,2]:=  Y[I]; AA[J,3]:=1.0;
        AA[J,4]:=+XX[I]; AA[J,5]:=+YY[I];
      end;
  NREF := J;    (* number of reference stars *)

  (* calculate plate constants a,b,c *)

  LSQFIT ( AA, NREF, 3, S );  A:=S[1]; B:=S[2]; C:=S[3];

  (* calculate plate constants d,e,f *)

  for I:=1 to NREF do  AA[I,4]:=AA[I,5]; (* copy column A[*,5]->A[*,4] *)
  LSQFIT ( AA, NREF, 3, S );  D:=S[1]; E:=S[2]; F:=S[3];

  (* calculate equatorial coordinates (and errors for reference stars) *)

  for I:=1 to NOBJ do
    begin
      XX[I] := A*X[I]+B*Y[I]+C;
      YY[I] := D*X[I]+E*Y[I]+F;
      STDEQU ( RA0,DEC0, XX[I],YY[I], RA_OBS,DEC_OBS );
      if NAME[I][1]='*' then (* find error in arcseconds *)
        begin
          D_RA  := (RA_OBS-RA[I])*CS(DEC[I]);
          D_DEC := (DEC_OBS-DEC[I]);
          DELTA[I] := 3600.0 * SQRT ( D_RA*D_RA + D_DEC*D_DEC );
        end;
      RA[I] := RA_OBS;  DEC[I] := DEC_OBS;
    end;

  (* focal length *)

  DET := A*E-D*B;
  FOC_LEN := 1.0/Sqrt(abs(DET));
  SCALE := ARC / FOC_LEN;

  (* output *)

  writeln (' Plate constants:' );
  writeln;
  writeln ('  a =',a:12:8,'  b =',b:12:8,'  c =',c:12:8);
  writeln ('  d =',d:12:8,'  e =',e:12:8,'  f =',f:12:8);
  writeln;
  writeln (' Effective focal length and image scale:');
  writeln;
  writeln ('  F =',FOC_LEN:9:2,' mm');
  writeln ('  m =',  SCALE:7:2,' "/mm');
  writeln;
  writeln (' Coordinates:');
  writeln;
  writeln (' Name':11, 'x':9,'y':7,'X':8,'Y':8,
           'RA':12, 'Dec':13, 'Error':9 );
  writeln ('mm':20, 'mm':7, ' ':23,
           'h  m  s  ', 'o  ''  " ':12, ' " ':6 );

  for I:=1 to NOBJ do
    begin
      WRITE('  '); for K:=1 to NAME_LENGTH do WRITE(NAME[I][K]);
      WRITE(X[I]:7:1,Y[I]:7:1,XX[I]:9:4,YY[I]:8:4);
      DMS(RA[I]/15.0,DEG,MIN,SEC); WRITE(DEG:5,MIN:3,SEC:6:2);
      DMS(DEC[I],DEG,MIN,SEC); WRITE(DEG:4,MIN:3,SEC:5:1);
      if NAME[I][1]='*' then WRITE(DELTA[I]:6:1);
      writeln;
    end;
  writeln;

end. (* FOTO *)

(*-----------------------------------------------------------------------*)
