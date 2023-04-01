program AstromifD;

uses
  Vcl.Forms,
  fAstro_d in 'source\interface\fAstro_d.pas' {frmAstroD},
  fAbout_d in 'source\interface\fAbout_d.pas' {frmAbout},
  uRandom in 'source\code\uRandom.pas',
  uSolarsys in 'source\code\uSolarsys.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmAstroD, frmAstroD);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
