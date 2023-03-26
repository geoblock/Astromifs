program AstromifD;

uses
  Vcl.Forms,
  fAstro_d in 'source\interface\fAstro_d.pas' {frmAstroD},
  Astro.Matlib in 'source\code\Astro.Matlib.pas',
  Astro.Sphlib in 'source\code\Astro.Sphlib.pas',
  Astro.Pnulib in 'source\code\Astro.Pnulib.pas',
  Astro.DiffEq in 'source\code\Astro.DiffEq.pas',
  Astro.Kepler in 'source\code\Astro.Kepler.pas',
  Astro.Moon in 'source\code\Astro.Moon.pas',
  Astro.Sunlib in 'source\code\Astro.Sunlib.pas',
  Astro.Planets in 'source\code\Astro.Planets.pas',
  Astro.Phylib in 'source\code\Astro.Phylib.pas',
  Astro.Timlib in 'source\code\Astro.Timlib.pas',
  fAbout_d in 'source\interface\fAbout_d.pas' {frmAbout};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmAstroD, frmAstroD);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
