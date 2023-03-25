program Astros;

uses
  Vcl.Forms,
  fAstros in 'source\vcl\fAstros.pas' {frmQosmos},
  Astro.Matlib in 'source\Astro.Matlib.pas',
  Astro.Sphlib in 'source\Astro.Sphlib.pas',
  Astro.Pnulib in 'source\Astro.Pnulib.pas',
  Astro.DiffEq in 'source\Astro.DiffEq.pas',
  Astro.Kepler in 'source\Astro.Kepler.pas',
  Astro.Moon in 'source\Astro.Moon.pas',
  Astro.Sunlib in 'source\Astro.Sunlib.pas',
  Astro.Planets in 'source\Astro.Planets.pas',
  Astro.Phylib in 'source\Astro.Phylib.pas',
  Astro.Timlib in 'source\Astro.Timlib.pas',
  fAbout in 'source\vcl\fAbout.pas' {frmAbout};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmQosmos, frmQosmos);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
