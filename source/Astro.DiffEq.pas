unit Astro.DiffEq;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Unit DELIB:                                                           *)
(*                                                                       *)
(*  Pascal implementation of the variable order variable stepsize        *)
(*  Adams-Bashforth-Moulton code DE for the numerical integration of     *)
(*  first order differential equations by Shampine and Gordon.           *)
(*                                                                       *)
(*  The original Fortran code is described in                            *)
(*                                                                       *)
(*  - Shampine L.F., Gordon M.K.: Computer solution of ordinary          *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*  - Shampine L.F., Gordon M.K.: Computer Loesung gewoehnlicher         *)
(*    Differentialgleichungen; Vieweg Verlag, Braunschweig (1984).       *)
(*  - Shampine L.F., Watts H.A.: DEPAC - Design of a user oriented       *)
(*    package of ODE solvers; Sand79-2374, Sandia Laboratories (1979).   *)
(*                                                                       *)
(*  The following functions and procedures are required to run DE:       *)
(*                                                                       *)
(*    Interp           Interpolation routine                             *)
(*    STEP            Basic integration routine                          *)
(*    DE              User called driver routine                         *)
(*                                                                       *)
(*  Furthermore, various type and constant definitions are required that *)
(*  have to be included into the declaration part of the main program.   *)
(*                                                                       *)
(*  The constant DE_EQN specifies the maximum dimension of the system of *)
(*  ordinary differential equations that can be integrated. TWOU and     *)
(*  FOURU have to be defined as two and four times the accuracy of the   *)
(*  Double data type.                                                    *)
(*                                                                       *)
(*    const DE_EQN = 6;  TWOU = 4.0E-12;  FOURU = 8.0E-12;               *)
(*                                                                       *)
(*  The DE_EQN_VECTOR type defines a vector of dimension DE_EQN for      *)
(*  state variables, while DE_12_VECTOR, DE_13_VECTOR and DE_PHI_VECTOR  *)
(*  are array types used to define various auxiliary variables.          *)
(*  DE_WORKSPACE_RECORD finally defines a record structure used to hold  *)
(*  all data that have to be stored globally between successive calls of *)
(*  the DE procedure.                                                    *)
(*                                                                       *)
(*    type  DE_EQN_VECTOR = array[1..DE_EQN] of Double;                  *)
(*          DE_12_VECTOR  = array[1..12] of Double;                      *)
(*          DE_13_VECTOR  = array[1..13] of Double;                      *)
(*          DE_PHI_ARRAY  = array[1..DE_EQN,1..16] of Double;            *)
(*          DE_WORKSPACE_RECORD = record                                 *)
(*                            YY,WT,P,YP,YPOUT     :  DE_EQN_VECTOR;     *)
(*                            PHI                  :  DE_PHI_ARRAY;      *)
(*                            ALPHA,BETA,V,W,PSI   :  DE_12_VECTOR;      *)
(*                            SIG,G                :  DE_13_VECTOR;      *)
(*                            X,H,HOLD,TOLD,DELSGN :  Double;            *)
(*                            NS,K,KOLD,ISNOLD     :  Integer;           *)
(*                            PHASE1,START,NORND   :  Boolean;           *)
(*                          end;                                         *)
(*                                                                       *)
(*                                                                       *)
(*  The function to be integrated must be declared as a procedure of     *)
(*  the form procedure ( var X: Double; var Y,YP: DE_EQN_VECTOR ); and   *)
(*  the same definition must be used in the procedure header of STEP     *)
(*  and DE.                                                              *)
(*                                                                       *)
(*  Using function EPSMACH the machine accuracy U=EPSMACH can be         *)
(*  computed to determine the proper values of constants TWOU and FOURU  *)
(*  for other compilers or floating point data types.                    *)
(*                                                                       *)
(*  Delphi adaption:                                                     *)
(*                                                                       *)
(*  Delphi uses special data types to pass a function or procedure       *)
(*  as formal argument of another function or procedure. For this        *)
(*  purpose the following procedure type is declared and used in the     *)
(*  specification of STEP and DE:                                        *)
(*                                                                       *)
(*      DE_FUNC = procedure ( X: Double; var Y,YP: DE_EQN_VECTOR );      *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)

interface

uses
  System.Math;

(* --------------------------------------------------------------------- *)
(* Sample dimensions and type declarations for use with DE,STEP,Interp *)
(* --------------------------------------------------------------------- *)

const
  DE_EQN = 6; (* Maximum number of ODEs to be integrated *)
  TWOU = 4.0E-12; (* Turbo Pascal Double data type accuracy *)
  FOURU = 8.0E-12; (* u=2**(-39)=1.81E-12 *)

type
  DE_EQN_VECTOR = array [1 .. DE_EQN] of Double;
  DE_FUNC = procedure(X: Double; Y: DE_EQN_VECTOR; var YP: DE_EQN_VECTOR);

  DE_12_VECTOR = array [1 .. 12] of Double;
  DE_13_VECTOR = array [1 .. 13] of Double;
  DE_PHI_ARRAY = array [1 .. DE_EQN, 1 .. 16] of Double;

  DE_WORKSPACE_RECORD = record
    YY, WT, P, YP, YPOUT: DE_EQN_VECTOR;
    PHI: DE_PHI_ARRAY;
    ALPHA, BETA, V, W, PSI: DE_12_VECTOR;
    SIG, G: DE_13_VECTOR;
    X, H, HOLD, TOLD, DELSGN: Double;
    NS, K, KOLD, ISNOLD: Integer;
    PHASE1, START, NORND: Boolean;
  end;

(* --------------------------------------------------------------------- *)
(*                          Procedures and functions                     *)
(* --------------------------------------------------------------------- *)

(*-----------------------------------------------------------------------*)
(* STEP                                                                  *)
(*                                                                       *)
(* Basic integration routine for the variable order and stepsize code    *)
(* of Shampine and Gordon.                                               *)
(*                                                                       *)
(*   X                                                                   *)
(*   Y                                                                   *)
(*   F       function to be integrated                                   *)
(*   NEQN    Number of differential equations                            *)
(*   H                                                                   *)
(*   EPS                                                                 *)
(*   WT                                                                  *)
(*   START                                                               *)
(*   HOLD                                                                *)
(*   K                                                                   *)
(*   KOLD                                                                *)
(*   CRASH   True if STEP was aborted (tolerances too small)             *)
(*   PHI                                                                 *)
(*   P                                                                   *)
(*   YP                                                                  *)
(*   PSI                                                                 *)
(*   ALPHA                                                               *)
(*   BETA                                                                *)
(*   SIG                                                                 *)
(*   V                                                                   *)
(*   W                                                                   *)
(*   G                                                                   *)
(*   PHASE1                                                              *)
(*   NS                                                                  *)
(*   NORND                                                               *)
(*                                                                       *)
(* Note:                                                                 *)
(*                                                                       *)
(* Replace "F:DE_FUNC;" by                                               *)
(*   "procedure F(X: double;Y:DE_EQN_VECTOR;var YP:DE_EQN_VECTOR);"         *)
(* for use with standard Pascal.                                         *)
(*                                                                       *)
(* A detailed description of the method is given in                      *)
(*                                                                       *)
(*    L. F. Shampine, M. K. Gordon; Computer Solution of ordinary        *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure STEP ( var X              :  Double;
                 var Y              :  DE_EQN_VECTOR;
                 F                  :  DE_FUNC;
                 var NEQN           :  Integer;
                 var H,EPS          :  Double;
                 var WT             :  DE_EQN_VECTOR;
                 var START          :  Boolean;
                 var HOLD           :  Double;
                 var K,KOLD         :  Integer;
                 var CRASH          :  Boolean;
                 var PHI            :  DE_PHI_ARRAY;
                 var P,YP           :  DE_EQN_VECTOR;
                 var PSI,ALPHA,BETA :  DE_12_VECTOR;
                 var SIG            :  DE_13_VECTOR;
                 var V,W            :  DE_12_VECTOR;
                 var G              :  DE_13_VECTOR;
                 var PHASE1         :  Boolean;
                 var NS             :  Integer;
                 var NORND          :  Boolean);


(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* Interp: Interpolation routine for use with DE and STEP                 *)
(*                                                                       *)
(*   X      Independent variable                                         *)
(*   Y      Solution at X                                                *)
(*   XOUT   Independent variable for which the solution is requested     *)
(*   YOUT   Interpolated solution at XOUT                                *)
(*   YPOUT  Derivative at XOUT                                           *)
(*   NEQN   Number of differential equations                             *)
(*   KOLD   Auxiliary quantity (most recently used order)                *)
(*   PHI    Auxiliary quantity                                           *)
(*   PSI    Auxiliary quantity                                           *)
(*                                                                       *)
(* STEP approximates the solution of the differential equation near X by *)
(* a polynomial, which is evaluated by Interp.                            *)
(*                                                                       *)
(* A detailed description of the method is given in                      *)
(*                                                                       *)
(*    L. F. Shampine, M. K. Gordon; Computer Solution of ordinary        *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*                                                                       *)
  (* ----------------------------------------------------------------------- *)
procedure Interp(X: Double; Y: DE_EQN_VECTOR; XOUT: Double; var YOUT, YPOUT: DE_EQN_VECTOR;
  NEQN, KOLD: Integer; PHI: DE_PHI_ARRAY; PSI: DE_12_VECTOR);

(*-----------------------------------------------------------------------*)
(* EPSMACH computes the machine accuracy u (1.0+u>1.0,1.0+u/2=1.0)       *)
(*-----------------------------------------------------------------------*)
function EPSMACH: Double;

(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* DE                                                                    *)
(*                                                                       *)
(* Driver routine for the variable order and step size Adams-Bashforth-  *)
(* Moulton method of Shampine and Gordon.                                *)
(*                                                                       *)
(*   F       function to be integrated                                   *)
(*   NEQN    Number of differential equations                            *)
(*   Y       State vector                                                *)
(*   T       Time                                                        *)
(*   TOUT    Output time                                                 *)
(*   RELERR  Relative tolerance                                          *)
(*   ABSERR  Absolute tolerance                                          *)
(*   IFLAG   Return code                                                 *)
(*   WORK    Work area for global storage between subsequent calls       *)
(*                                                                       *)
(* Note:                                                                 *)
(*                                                                       *)
(* Subroutine DE is used for the numerical integration of an ordinary    *)
(* differential equation y'=f(x,y) of first order, which must be         *)
(* provided as a subroutine of the form                                  *)
(*   procedure F ( X: Double; Y: DE_EQN_VECTOR; var YP: DE_EQN_VECTOR );   *)
(* Aside from DE itself the subroutines STEP and Interp as well as some   *)
(* global constant and type definitions are required, which are          *)
(* described in detail in DELIB.PAS.                                     *)
(*                                                                       *)
(* On the first call of DE the initial time is passed in T, while Y      *)
(* contains the corresponding state vector. Furthermore the desired      *)
(* output time TOUT and the required relative and absolute accuracy      *)
(* RELERR and ABSERR need to be specified. The status variable IFLAG     *)
(* is set to 1 to inform DE that the integration of a new problem is     *)
(* started.                                                              *)
(*                                                                       *)
(* In normal cases the variable Y subsequently contains the state vector *)
(* at time TOUT. Simultaneously T is assigned the value TOUT, while the  *)
(* status flag IFLAg is set to 2 (successful step). In all continuation  *)
(* steps it is only required for the user to specify a new value TOUT,   *)
(* while all remaining variables (most notably IFLAG) are left           *)
(* unchanged. IFLAG needs only be reset to one in case of a new problem  *)
(* or when changing the direction of integration.                        *)
(*                                                                       *)
(* If DE returns a value of IFLAG different from 2, one of the following *)
(* cases may have occured:                                               *)
(*                                                                       *)
(*   IFLAG = 3  TOUT has not been reached, since the requested           *)
(*              tolerance RELERR or ABSERR was too small. Both           *)
(*              values have been increased for subsequent                *)
(*              computations.                                            *)
(*   IFLAG = 4  TOUT has not been reached, since more than               *)
(*              MAXNUM=500 steps have been required internally.          *)
(*   IFLAG = 5  TOUT has not been reached, since the differential        *)
(*              equation appears to be stiff. This should not happen     *)
(*              when applying DE to problems in celestial mechanics.     *)
(*   IFLAG = 6  Illegal input parameter (e.g. T=TOUT)                    *)
(*                                                                       *)
(* In all cases except IFLAG=6 the integration may simply be continued   *)
(* by calling DE again without changing any parmeter. The interruption   *)
(* is mainly intended to call the user's attention to possible problems  *)
(* that may have occured and e.g. avoid infinite integrations.           *)
(*                                                                       *)
(* The work array WORK is used to store various quantities between       *)
(* subsequent calls of DE and is not changed by the user.                *)
(*                                                                       *)
(* Note that "F:DE_FUNC;" in the specification of DE has to be replaced  *)
(* by "procedure F(X: double;Y:DE_EQN_VECTOR;var YP:DE_EQN_VECTOR);" for    *)
(* use with standard Pascal.                                             *)
(*                                                                       *)
(* A detailed description of the method and its Fortran implementation   *)
(* is given in                                                           *)
(*                                                                       *)
(*  - Shampine L.F., Gordon M.K.: Computer solution of ordinary          *)
(*    Differential Equations; Freeman and Comp., San Francisco (1975).   *)
(*  - Shampine L.F., Watts H.A.: DEPAC - Design of a user oriented       *)
(*    package of ODE solvers; Sand79-2374, Sandia Laboratories (1979).   *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure DE(F: DE_FUNC; NEQN: Integer; var Y: DE_EQN_VECTOR; var T, TOUT: Double;
  var RELERR, ABSERR: Double; var IFLAG: Integer; var WORK: DE_WORKSPACE_RECORD);


//==============================================================
implementation
//==============================================================

(*-----------------------------------------------------------------------*)

procedure STEP(var X: Double; var Y: DE_EQN_VECTOR; F: DE_FUNC; var NEQN: Integer;
  var H, EPS: Double; var WT: DE_EQN_VECTOR; var START: Boolean; var HOLD: Double;
  var K, KOLD: Integer; var CRASH: Boolean; var PHI: DE_PHI_ARRAY; var P, YP: DE_EQN_VECTOR;
  var PSI, ALPHA, BETA: DE_12_VECTOR; var SIG: DE_13_VECTOR; var V, W: DE_12_VECTOR;
  var G: DE_13_VECTOR; var PHASE1: Boolean; var NS: Integer; var NORND: Boolean);

label 999; (* Error exit at end of procedure *)

var
  I, IFAIL, IM1, IP1, IQ, J, KM1, KM2, KNEW, KP1, KP2: Integer;
  L, LIMIT1, LIMIT2, NSM2, NSP1, NSP2: Integer;
  ABSH, ERK, ERKM1, ERKM2, ERKP1, ERR, HNEW: Double;
  P5EPS, R, REALI, REALNS, RHO, ROUND, SUM, TAU, TEMP1: Double;
  TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, XOLD: Double;
  TWO, GSTR: DE_13_VECTOR;
  SUCCESS: Boolean;

  (* SIGN: computes sign(B)*abs(A) *)
  (*sub*)function SIGN(A, B: Double): Double;
  begin
    if (B >= 0.0) then
      SIGN := abs(A)
    else
      SIGN := -abs(A);
  end;

  (* MIN: computes the minimum of A and B *)
  (*sub*)function MIN(A, B: Double): Double;
  begin
    if A < B then
      MIN := A
    else
      MIN := B;
  end;

  (* MAX: computes the maximum of A and B *)
  (*sub*)function MAX(A, B: Double): Double;
  begin
    if A > B then
      MAX := A
    else
      MAX := B;
  end;

begin (* STEP *)

  (* explicit initialization of arrays TWO and GSTR *)

  TWO[1] := 2.0;
  for I := 2 to 13 do
    TWO[I] := 2.0 * TWO[I - 1];

  GSTR[1] := 0.5;
  GSTR[2] := 0.0833;
  GSTR[3] := 0.0417;
  GSTR[4] := 0.0264;
  GSTR[5] := 0.0188;
  GSTR[6] := 0.0143;
  GSTR[7] := 0.0114;
  GSTR[8] := 0.00936;
  GSTR[9] := 0.00789;
  GSTR[10] := 0.00679;
  GSTR[11] := 0.00592;
  GSTR[12] := 0.00524;
  GSTR[13] := 0.00468;

  (* *)
  (* Begin block 0 *)
  (* *)
  (* Check if step size or error tolerance is too small for machine *)
  (* precision.  If first step, initialize PHI array and estimate a *)
  (* starting step size. If step size is too small, determine an *)
  (* acceptable one. *)
  (* *)

  if (abs(H) < FOURU * abs(X)) then
  begin
    H := SIGN(FOURU * abs(X), H);
    CRASH := True;
    goto 999; (* Exit *)
  end;

  P5EPS := 0.5 * EPS;
  CRASH := False;
  G[1] := 1.0;
  G[2] := 0.5;
  SIG[1] := 1.0;

  IFAIL := 0;

  (* If error tolerance is too small, increase it to an *)
  (* acceptable value. *)

  ROUND := 0.0;
  for L := 1 to NEQN do
    ROUND := ROUND + (Y[L] * Y[L]) / (WT[L] * WT[L]);
  ROUND := TWOU * Sqrt(ROUND);
  if P5EPS < ROUND then
  begin
    EPS := 2.0 * ROUND * (1.0 + FOURU);
    CRASH := True;
    goto 999;
  end;

  if (START) then
  begin
    (* Initialize. Compute appropriate step size for first step. *)
    F(X, Y, YP);
    SUM := 0.0;
    for L := 1 to NEQN do
    begin
      PHI[L, 1] := YP[L];
      PHI[L, 2] := 0.0;
      SUM := SUM + (YP[L] * YP[L]) / (WT[L] * WT[L]);
    end;
    SUM := Sqrt(SUM);
    ABSH := abs(H);
    if (EPS < 16.0 * SUM * H * H) then
      ABSH := 0.25 * Sqrt(EPS / SUM);
    H := SIGN(MAX(ABSH, FOURU * abs(X)), H);
    HOLD := 0.0;
    K := 1;
    KOLD := 0;
    START := False;
    PHASE1 := True;
    NORND := True;
    if (P5EPS <= 100.0 * ROUND) then
    begin
      NORND := False;
      for L := 1 to NEQN do
        PHI[L, 15] := 0.0;
    end;
  end;

  (* *)
  (* End block 0 *)
  (* *)

  (* *)
  (* Repeat blocks 1, 2 (and 3) until step is successful *)
  (* *)

  repeat
    (* *)
    (* Begin block 1 *)
    (* *)
    (* Compute coefficients of formulas for this step. Avoid computing *)
    (* those quantities not changed when step size is not changed. *)
    (* *)
    KP1 := K + 1;
    KP2 := K + 2;
    KM1 := K - 1;
    KM2 := K - 2;

    (* NS is the number of steps taken with size H, including the *)
    (* current one.  When K<NS, no coefficients change. *)
    if (H <> HOLD) then
      NS := 0;
    if (NS <= KOLD) then
      NS := NS + 1;
    NSP1 := NS + 1;

    if (K >= NS) then
    begin
      (* Compute those components of ALPHA[*],BETA[*],PSI[*],SIG[*] *)
      (* which are changed *)
      BETA[NS] := 1.0;
      REALNS := NS;
      ALPHA[NS] := 1.0 / REALNS;
      TEMP1 := H * REALNS;
      SIG[NSP1] := 1.0;
      if (K >= NSP1) then
        for I := NSP1 to K do
        begin
          IM1 := I - 1;
          TEMP2 := PSI[IM1];
          PSI[IM1] := TEMP1;
          BETA[I] := BETA[IM1] * PSI[IM1] / TEMP2;
          TEMP1 := TEMP2 + H;
          ALPHA[I] := H / TEMP1;
          REALI := I;
          SIG[I + 1] := REALI * ALPHA[I] * SIG[I];
        end;
      PSI[K] := TEMP1;

      (* Compute coefficients G[*]; initialize V[*] and set W[*]. *)

      if (NS > 1) then
      begin
        (* If order was raised, update diagonal part of V[*] *)
        if (K > KOLD) then
        begin
          TEMP4 := K * KP1;
          V[K] := 1.0 / TEMP4;
          NSM2 := NS - 2;
          for J := 1 to NSM2 do
          begin
            I := K - J;
            V[I] := V[I] - ALPHA[J + 1] * V[I + 1];
          end;
        end;

        (* Update V[*] and set W[*] *)
        LIMIT1 := KP1 - NS;
        TEMP5 := ALPHA[NS];
        for IQ := 1 to LIMIT1 do
        begin
          V[IQ] := V[IQ] - TEMP5 * V[IQ + 1];
          W[IQ] := V[IQ];
        end;
        G[NSP1] := W[1];

      end
      else
        for IQ := 1 to K do
        begin
          TEMP3 := IQ * (IQ + 1);
          V[IQ] := 1.0 / TEMP3;
          W[IQ] := V[IQ];
        end;

      (* Compute the G[*] in the work vector W[*] *)
      NSP2 := NS + 2;
      if (KP1 >= NSP2) then
        for I := NSP2 to KP1 do
        begin
          LIMIT2 := KP2 - I;
          TEMP6 := ALPHA[I - 1];
          for IQ := 1 to LIMIT2 do
            W[IQ] := W[IQ] - TEMP6 * W[IQ + 1];
          G[I] := W[1];
        end;
    end; (* if K>=NS *)

    (* *)
    (* End block 1 *)
    (* *)

    (* *)
    (* Begin block 2 *)
    (* *)
    (* Predict a solution P[*], evaluate derivatives using predicted *)
    (* solution, estimate local error at order K and errors at orders *)
    (* K, K-1, K-2 as if constant step size were used. *)
    (* *)

    (* Change PHI to PHI star *)
    if (K >= NSP1) then
      for I := NSP1 to K do
      begin
        TEMP1 := BETA[I];
        for L := 1 to NEQN do
          PHI[L, I] := TEMP1 * PHI[L, I];
      end;

    (* Predict solution and differences *)
    for L := 1 to NEQN do
    begin
      PHI[L, KP2] := PHI[L, KP1];
      PHI[L, KP1] := 0.0;
      P[L] := 0.0;
    end;
    for J := 1 to K do
    begin
      I := KP1 - J;
      IP1 := I + 1;
      TEMP2 := G[I];
      for L := 1 to NEQN do
      begin
        P[L] := P[L] + TEMP2 * PHI[L, I];
        PHI[L, I] := PHI[L, I] + PHI[L, IP1];
      end;
    end;
    if (NORND) then
      for L := 1 to NEQN do
        P[L] := Y[L] + H * P[L]
    else
      for L := 1 to NEQN do
      begin
        TAU := H * P[L] - PHI[L, 15];
        P[L] := Y[L] + TAU;
        PHI[L, 16] := (P[L] - Y[L]) - TAU;
      end;
    XOLD := X;
    X := X + H;
    ABSH := abs(H);
    F(X, P, YP);

    (* Estimate errors at orders K,K-1,K-2 *)
    ERKM2 := 0.0;
    ERKM1 := 0.0;
    ERK := 0.0;
    for L := 1 to NEQN do
    begin
      TEMP3 := 1.0 / WT[L];
      TEMP4 := YP[L] - PHI[L, 1];
      if (KM2 > 0) then
        ERKM2 := ERKM2 + ((PHI[L, KM1] + TEMP4) * TEMP3) * ((PHI[L, KM1] + TEMP4) * TEMP3);
      if (KM2 >= 0) then
        ERKM1 := ERKM1 + ((PHI[L, K] + TEMP4) * TEMP3) * ((PHI[L, K] + TEMP4) * TEMP3);
      ERK := ERK + (TEMP4 * TEMP3) * (TEMP4 * TEMP3);
    end;
    if (KM2 > 0) then
      ERKM2 := ABSH * SIG[KM1] * GSTR[KM2] * Sqrt(ERKM2);
    if (KM2 >= 0) then
      ERKM1 := ABSH * SIG[K] * GSTR[KM1] * Sqrt(ERKM1);
    TEMP5 := ABSH * Sqrt(ERK);
    ERR := TEMP5 * (G[K] - G[KP1]);
    ERK := TEMP5 * SIG[KP1] * GSTR[K];
    KNEW := K;

    (* Test if order should be lowered *)
    if (KM2 > 0) then
      if (MAX(ERKM1, ERKM2) <= ERK) then
        KNEW := KM1;
    if (KM2 = 0) then
      if (ERKM1 <= 0.5 * ERK) then
        KNEW := KM1;

    (* *)
    (* End block 2 *)
    (* *)

    (* *)
    (* If step is successful continue with block 4, otherwise repeat *)
    (* blocks 1 and 2 after executing block 3 *)
    (* *)

    SUCCESS := (ERR <= EPS);
    if (NOT SUCCESS) then
    begin
      (* *)
      (* Begin block 3 *)
      (* *)
      (* The step is unsuccessful. Restore X, PHI[*,*], PSI[*]. If *)
      (* 3rd consecutive failure, set order to 1. If step fails more *)
      (* than 3 times, consider an optimal step size. Double error *)
      (* tolerance and return if estimated step size is too small *)
      (* for machine precision. *)
      (* *)

      (* Restore X, PHI[*,*] and PSI[*] *)
      PHASE1 := False;
      X := XOLD;
      for I := 1 to K do
      begin
        TEMP1 := 1.0 / BETA[I];
        IP1 := I + 1;
        for L := 1 to NEQN do
          PHI[L, I] := TEMP1 * (PHI[L, I] - PHI[L, IP1]);
      end;
      if (K >= 2) then
        for I := 2 to K do
          PSI[I - 1] := PSI[I] - H;

      (* On third failure, set order to one. *)
      (* Thereafter, use optimal step size *)
      IFAIL := IFAIL + 1;
      TEMP2 := 0.5;
      if (IFAIL > 3) then
        if (P5EPS < 0.25 * ERK) then
          TEMP2 := Sqrt(P5EPS / ERK);
      if (IFAIL >= 3) then
        KNEW := 1;
      H := TEMP2 * H;
      K := KNEW;
      if (abs(H) < FOURU * abs(X)) then
      begin
        CRASH := True;
        H := SIGN(FOURU * abs(X), H);
        EPS := EPS + EPS;
        goto 999; (* Exit *)
      end;

      (* *)
      (* End block 3, return to start of block 1 *)
      (* *)

    end; (* end if(successful) *)
  until (SUCCESS);

  (* *)
  (* Begin block 4 *)
  (* *)
  (* The step is successful. Correct the predicted solution, evaluate *)
  (* the derivatives using the corrected solution and update the *)
  (* differences. Determine best order and step size for next step. *)
  (* *)

  KOLD := K;
  HOLD := H;

  (* Correct and evaluate *)
  TEMP1 := H * G[KP1];
  if (NORND) then
    for L := 1 to NEQN do
      Y[L] := P[L] + TEMP1 * (YP[L] - PHI[L, 1])
  else
    for L := 1 to NEQN do
    begin
      RHO := TEMP1 * (YP[L] - PHI[L, 1]) - PHI[L, 16];
      Y[L] := P[L] + RHO;
      PHI[L, 15] := (Y[L] - P[L]) - RHO;
    end;
  F(X, Y, YP);

  (* Update differences for next step *)
  for L := 1 to NEQN do
  begin
    PHI[L, KP1] := YP[L] - PHI[L, 1];
    PHI[L, KP2] := PHI[L, KP1] - PHI[L, KP2];
  end;
  for I := 1 to K do
    for L := 1 to NEQN do
      PHI[L, I] := PHI[L, I] + PHI[L, KP1];

  (* Estimate error at order K+1 unless *)
  (* - in first phase when always raise order, *)
  (* - already decided to lower order, *)
  (* - step size not constant so estimate unreliable *)

  ERKP1 := 0.0;
  if ((KNEW = KM1) or (K = 12)) then
    PHASE1 := False;

  if (PHASE1)
  then
  begin
    K := KP1;
    ERK := ERKP1;
  end
  else
  begin
    if (KNEW = KM1)
    then (* lower order *)
    begin
      K := KM1;
      ERK := ERKM1;
    end
    else
      if (KP1 <= NS) then
    begin
      for L := 1 to NEQN do
        ERKP1 := ERKP1 + (PHI[L, KP2] / WT[L]) * (PHI[L, KP2] / WT[L]);
      ERKP1 := ABSH * GSTR[KP1] * Sqrt(ERKP1);
      (* Using estimated error at order K+1, determine *)
      (* appropriate order for next step *)
      if (K > 1) then
        if (ERKM1 <= MIN(ERK, ERKP1)) then (* lower order *)
        begin
          K := KM1;
          ERK := ERKM1
        end
        else
        begin
          if ((ERKP1 < ERK) and (K <> 12)) then
          (* raise order *)
          begin
            K := KP1;
            ERK := ERKP1;
          end;
        end
      else if (ERKP1 < 0.5 * ERK) then (* raise order *)
      begin
        (* Here ERKP1 < ERK < max(ERKM1,ERKM2) else *)
        (* order would have been lowered in block 2. *)
        (* Thus order is to be raised *)
        K := KP1;
        ERK := ERKP1;
      end;
      (* end if  K>1 *)
    end; (* if KP1<=NS *)

    (* end if KNEW=KM1 *)

  end; (* if PHASE1 *)

  (* With new order determine appropriate step size for next step *)
  if ((PHASE1) or (P5EPS >= ERK * TWO[K + 1])) then
    HNEW := H + H
  else
  begin
    if (P5EPS < ERK) then
    begin
      TEMP2 := K + 1;
      R := EXP(LN(P5EPS / ERK) * (1.0 / TEMP2));
      HNEW := ABSH * MAX(0.5, MIN(0.9, R));
      HNEW := SIGN(MAX(HNEW, FOURU * abs(X)), H);
    end
    else
      HNEW := H;
  end;

  H := HNEW;
  (* *)
  (* End block 4 *)
  (* *)

999:
end; (* STEP *)


(* ----------------------------------------------------------------------- *)
procedure Interp(X: Double; Y: DE_EQN_VECTOR; XOUT: Double; var YOUT, YPOUT: DE_EQN_VECTOR;
  NEQN, KOLD: Integer; PHI: DE_PHI_ARRAY; PSI: DE_12_VECTOR);
var
  I, J, JM1, KI, KIP1, L, LIMIT1: Integer;
  ETA, GAMMA, HI, PSIJM1: Double;
  TEMP1, TEMP2, TEMP3, TERM: Double;
  G, RHO, W: DE_13_VECTOR;

begin
  G[1] := 1.0;
  RHO[1] := 1.0;

  HI := XOUT - X;
  KI := KOLD + 1;
  KIP1 := KI + 1;

  (* Initialize W[*] for computing G[*] *)
  for I := 1 to KI do
  begin
    TEMP1 := I;
    W[I] := 1.0 / TEMP1;
  end;
  TERM := 0.0;

  (* compute G[*] *)
  for J := 2 to KI do
  begin
    JM1 := J - 1;
    PSIJM1 := PSI[JM1];
    GAMMA := (HI + TERM) / PSIJM1;
    ETA := HI / PSIJM1;
    LIMIT1 := KIP1 - J;
    for I := 1 to LIMIT1 do
      W[I] := GAMMA * W[I] - ETA * W[I + 1];
    G[J] := W[1];
    RHO[J] := GAMMA * RHO[JM1];
    TERM := PSIJM1;
  end;

  (* Interpolate for the solution YOUT and for *)
  (* the derivative of the solution YPOUT *)
  for L := 1 to NEQN do
  begin
    YPOUT[L] := 0.0;
    YOUT[L] := 0.0;
  end;

  for J := 1 to KI do
  begin
    I := KIP1 - J;
    TEMP2 := G[I];
    TEMP3 := RHO[I];
    for L := 1 to NEQN do
    begin
      YOUT[L] := YOUT[L] + TEMP2 * PHI[L, I];
      YPOUT[L] := YPOUT[L] + TEMP3 * PHI[L, I];
    end;
  end;

  for L := 1 to NEQN do
    YOUT[L] := Y[L] + HI * YOUT[L];
end; (* Interp *)

//----------------------------------------------------------------------------

procedure DE(F: DE_FUNC; (* Turbo Pascal *)
  NEQN: Integer; var Y: DE_EQN_VECTOR; var T, TOUT: Double; var RELERR, ABSERR: Double;
  var IFLAG: Integer; var WORK: DE_WORKSPACE_RECORD);

label 99; (* Error exit at end of procedure *)
const
  MAXNUM = 500; (* Maximum number of steps *)

var
  STIFF, CRASH: Boolean;
  ISN, NOSTEP, KLE4, L: Integer;
  RELEPS, ABSEPS, TEND, ABSDEL, DEL, EPS: Double;

  (* SIGN: computes sign(B)*abs(A) *)
  (*sub*)function SIGN(A, B: Double): Double;
  begin
    if (B >= 0.0) then
      SIGN := abs(A)
    else
      SIGN := -abs(A);
  end;

  (* MIN: computes the minimum of A and B *)
  (*sub*)function MIN(A, B: Double): Double;
  begin
    if A < B then
      MIN := A
    else
      MIN := B;
  end;

(* MAX: computes the maximum of A and B *)

  (*sub*)function MAX(A, B: Double): Double;
  begin
    if A > B then
      MAX := A
    else
      MAX := B;
  end;

begin
  with WORK do (* Use short form of record component name *)
  begin
    (* Test for invalid parameters *)
    EPS := MAX(RELERR, ABSERR);
    ISN := Trunc(SIGN(1.1, IFLAG)); (* should be +/- 1 *)
    IFLAG := abs(IFLAG);
    if ((NEQN < 1) or (T = TOUT) or (RELERR < 0.0) or (ABSERR < 0.0) or (EPS = 0.0) or
      (IFLAG = 0) or (IFLAG > 5) or ((IFLAG <> 1) and (T <> TOLD))) then
    begin
      IFLAG := 6; (* Set error code *)
      goto 99; (* Exit *)
    end;

    (* Set interval limits and step counter *)
    DEL := TOUT - T;
    ABSDEL := abs(DEL);
    TEND := T + 10.0 * DEL;
    if (ISN < 0) then
      TEND := TOUT;
    NOSTEP := 0;
    KLE4 := 0;
    STIFF := False;
    RELEPS := RELERR / EPS;
    ABSEPS := ABSERR / EPS;

    if ((IFLAG = 1) or (ISNOLD < 0) or (DELSGN * DEL <= 0.0)) then
    begin
      (* Set independent and dependent variables X and YY[*] for *)
      (* steps. Set sign of integration direction. Initialize the *)
      (* step size. *)

      START := True;
      X := T;
      for L := 1 to NEQN do
        YY[L] := Y[L];
      DELSGN := SIGN(1.0, DEL);
      H := SIGN(MAX(FOURU * abs(X), abs(TOUT - X)), TOUT - X);
    end;
    repeat
      (* If already past output point, interpolate and return *)
      if (abs(X - T) >= ABSDEL) then
      begin
        Interp(X, YY, TOUT, Y, YPOUT, NEQN, KOLD, PHI, PSI);
        IFLAG := 2;
        T := TOUT;
        TOLD := T;
        ISNOLD := ISN;
        goto 99; (* Exit *)
      end;

      (* If cannot go past TSTOP and sufficiently close, *)
      (* extrapolate and return *)
      if ((ISN <= 0) and (abs(TOUT - X) < FOURU * abs(X))) then
      begin
        H := TOUT - X;
        F(X, YY, YP);
        for L := 1 to NEQN do
          Y[L] := YY[L] + H * YP[L];
        IFLAG := 2;
        T := TOUT;
        TOLD := T;
        ISNOLD := ISN;
        goto 99; (* Exit *)
      end;

      (* Monitor number of steps attempted *)

      if (NOSTEP >= MAXNUM) then
      (* a large amount of work has been expended *)
      begin
        IFLAG := ISN * 4;
        if (STIFF) then
          IFLAG := ISN * 5;
        for L := 1 to NEQN do
          Y[L] := YY[L];
        T := X;
        TOLD := T;
        ISNOLD := 1;
        goto 99; (* Exit *)
      end;

      (* Limit step size, set weight vector and take a step *)

      H := SIGN(MIN(abs(H), abs(TEND - X)), H);
      for L := 1 to NEQN do
        WT[L] := RELEPS * abs(YY[L]) + ABSEPS;

      STEP(X, YY, F, NEQN, H, EPS, WT, START, HOLD, K, KOLD, CRASH, PHI, P, YP, PSI, ALPHA, BETA,
        SIG, V, W, G, PHASE1, NS, NORND);

      if (CRASH) then (* Tolerances too small *)
      begin
        IFLAG := ISN * 3;
        RELERR := EPS * RELEPS;
        ABSERR := EPS * ABSEPS;
        for L := 1 to NEQN do
          Y[L] := YY[L];
        T := X;
        TOLD := T;
        ISNOLD := 1;
        goto 99; (* Exit *)
      end;

      (* Stiffness test: *)
      (* count number of consecutive steps taken with the *)
      (* order of the method being less or equal to four *)
      NOSTEP := NOSTEP + 1;
      KLE4 := KLE4 + 1;
      if (KOLD > 4) then
        KLE4 := 0;
      if (KLE4 >= 50) then
        STIFF := True;

    until False;
  end; (* with WORK do *)
99:
end;

(*-----------------------------------------------------------------------*)

function EPSMACH: Double;
var
  ONE, TWO, U: Double;
begin
  ONE := 1.0;
  TWO := 2.0;
  U := 1.0;
  repeat
    U := U / TWO;
  until ((ONE + U) = ONE);
  EPSMACH := TWO * U;
end;


end.
