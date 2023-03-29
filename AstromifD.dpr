program AstromifD;

uses
  Vcl.Forms,
  fAstro_d in 'source\interface\fAstro_d.pas' {frmAstroD},
  fAbout_d in 'source\interface\fAbout_d.pas' {frmAbout},
  astronomy in 'source\astronomy\astronomy.pas',
  Qod.Randomam in 'source\code\Qod.Randomam.pas',
  Qod.SolarSys in 'source\code\Qod.SolarSys.pas',
  Apc.DiffEq in 'source\apc\Apc.DiffEq.pas',
  Apc.Mathem in 'source\apc\Apc.Mathem.pas',
  Apc.Moon in 'source\apc\Apc.Moon.pas',
  Apc.Planets in 'source\apc\Apc.Planets.pas',
  Apc.Sun in 'source\apc\Apc.Sun.pas',
  Apc.Time in 'source\apc\Apc.Time.pas',
  Apc.Physic in 'source\apc\Apc.Physic.pas',
  Apc.Kepler in 'source\apc\Apc.Kepler.pas',
  Apc.Spheric in 'source\apc\Apc.Spheric.pas',
  Apc.PrecNut in 'source\apc\Apc.PrecNut.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmAstroD, frmAstroD);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
