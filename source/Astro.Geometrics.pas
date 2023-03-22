unit Astro.Geometrics;

interface

uses
  Winapi.OpenGL;

procedure SpheraGL(r: Real; gm1, gm2: Word);

(*
SpheraGL(0.9,720,360);   //1 sphere
SpheraGL(0.9, 90,720);   //1/2 sphere without vertices
SpheraGL(0.9,180,360):   //1/4 sphere
SpheraGL(0.9,360,360);   //1/2 sphere
SpheraGL(0.9,540,360);   //3/4 sphere
*)


implementation


procedure SpheraGL(r: Real; gm1, gm2: Word);
var
  g1, g2: Integer;
  xx, yy, zz, gr1, gr2: Real;
const
  p = pi / 360;
begin
  glPointSize(3);
  glBegin(GL_POINTS);
  for g1 := 0 to gm1 do
    for g2 := 0 to gm2 do
    begin
      gr1 := g1 * p;
      gr2 := g2 * p;
      xx := r * cos(gr1) * cos(gr2);
      yy := r * cos(gr1) * sin(gr2);
      zz := r * sin(gr1);
      glColor(xx, yy, zz);
      glVertex3f(xx, yy, zz);
    end;
  glEnd;
end;

end.
