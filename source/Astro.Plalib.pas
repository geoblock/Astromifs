unit Astro.Plalib;

interface

uses
  Astro.Matlib,
  Astro.Sphlib,
  Astro.Keplib;

type
  PlanetType = (Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto);


(*-----------------------------------------------------------------------*)
(*                                                                       *)
(* POSITION:                                                             *)
(*                                                                       *)
(*   Computes the position of a planet assuming Keplerian orbits. Mean   *)
(*   elements at epoch J2000 are used for Mercury to Mars; osculating    *)
(*   elements at epoch 1995/10/10 (JD 2450000.5) are used for Jupiter to *)
(*   Pluto. Relative accuracy approx. 0.001 between 1990 and 2000.       *)
(*                                                                       *)
(*   Planet  Name of the planet                                          *)
(*   T       Time in Julian centuries since J2000                        *)
(*   X,Y,Z   Ecliptic coordinates (equinox J2000)                        *)
(*                                                                       *)
(*-----------------------------------------------------------------------*)
procedure Position(Planet: PlanetType; T: Double; var X, Y, Z: Double);

//=================================================================
implementation
//=================================================================

procedure Position(Planet: PlanetType; T: Double; var X, Y, Z: Double);
var
  A, E, M, O, I, W, N: Double;
  XX, YY, VVX, VVY: Double;
  T0: Double;
  PQR: Double33;
begin
  (* Heliocentric ecliptic orbital elements for equinox J2000 *)
  case Planet of
    Mercury:
      begin
        A := 0.387099;
        E := 0.205634;
        M := 174.7947;
        N := 149472.6738;
        O := 48.331;
        I := 7.0048;
        W := 77.4552;
      end;
    Venus:
      begin
        A := 0.723332;
        E := 0.006773;
        M := 50.4071;
        N := 58517.8149;
        O := 76.680;
        I := 3.3946;
        W := 131.5718;
      end;
    Earth:
      begin
        A := 1.000000;
        E := 0.016709;
        M := 357.5256;
        N := 35999.3720;
        O := 174.876;
        I := 0.0000;
        W := 102.9400;
      end;
    Mars:
      begin
        A := 1.523692;
        E := 0.093405;
        M := 19.3879;
        N := 19140.3023;
        O := 49.557;
        I := 1.8496;
        W := 336.0590;
      end;
    Jupiter:
      begin
        A := 5.202437;
        E := 0.048402;
        M := 250.3274;
        N := 3035.2275;
        O := 100.4683;
        I := 1.3047;
        W := 15.7192;
      end;
    Saturn:
      begin
        A := 9.551712;
        E := 0.052340;
        M := 267.2465;
        N := 1219.6465;
        O := 113.6439;
        I := 2.4855;
        W := 90.9682;
      end;
    Uranus:
      begin
        A := 19.293108;
        E := 0.044846;
        M := 118.4320;
        N := 424.8150;
        O := 74.0903;
        I := 0.7733;
        W := 176.6152;
      end;
    Neptune:
      begin
        A := 30.257162;
        E := 0.007985;
        M := 292.4716;
        N := 216.3047;
        O := 131.7750;
        I := 1.7700;
        W := 3.0962;
      end;
    Pluto:
      begin
        A := 39.783607;
        E := 0.254351;
        M := 8.2304;
        N := 143.4629;
        O := 110.3865;
        I := 17.1201;
        W := 224.7424;
      end;
  end;
  case Planet of
    Mercury, Venus, Earth, Mars:
      T0 := 0.0;
    Jupiter, Saturn, Uranus, Neptune, Pluto:
      T0 := -0.042286105;
  end;
  M := M + N * (T - T0);

  (* Cartesian coordinates mean ecliptic and equinox J2000 *)
  GAUSVEC(O, I, W - O, PQR);
  ELLIP(M, A, E, XX, YY, VVX, VVY);
  ORBECL(XX, YY, PQR, X, Y, Z);
end; // Position

end.
