program AstromifD;

uses
  Vcl.Forms,
  fAstro_d in 'source\interface\fAstro_d.pas' {frmAstroD},
  fAbout_d in 'source\interface\fAbout_d.pas' {frmAbout},
  astronomy in 'source\astronomy\astronomy.pas',
  Qod.Randomam in 'source\code\Qod.Randomam.pas',
  Qod.SolarSys in 'source\code\Qod.SolarSys.pas',
  Apc.DiffEq in 'source\apc\Apc.DiffEq.pas',
  Apc.Kepler in 'source\apc\Apc.Kepler.pas',
  Apc.Matlib in 'source\apc\Apc.Matlib.pas',
  Apc.Moon in 'source\apc\Apc.Moon.pas',
  Apc.Phylib in 'source\apc\Apc.Phylib.pas',
  Apc.Planets in 'source\apc\Apc.Planets.pas',
  Apc.Pnulib in 'source\apc\Apc.Pnulib.pas',
  Apc.Sphlib in 'source\apc\Apc.Sphlib.pas',
  Apc.Sunlib in 'source\apc\Apc.Sunlib.pas',
  Apc.Timlib in 'source\apc\Apc.Timlib.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmAstroD, frmAstroD);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
