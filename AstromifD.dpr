program AstromifD;

uses
  Vcl.Forms,
  fAstro_d in 'source\interface\fAstro_d.pas' {frmAstroD},
  fAbout_d in 'source\interface\fAbout_d.pas' {frmAbout},
  Qod.Randomam in 'source\code\Qod.Randomam.pas',
  Qod.SolarSys in 'source\code\Qod.SolarSys.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmAstroD, frmAstroD);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
