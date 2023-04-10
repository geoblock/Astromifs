program AstromifD;

uses
  Vcl.Forms,
  fSonofon in 'source\interface\fSonofon.pas' {frmSonofon},
  fInitials in 'source\interface\fInitials.pas' {FormInitial},
  fOptions in 'source\interface\fOptions.pas' {FormOptions},
  fDialogs in 'source\interface\fDialogs.pas' {FormDialog},
  fAbout in 'source\interface\fAbout.pas' {frmAbout},
  fAstromifs in 'source\interface\fAstromifs.pas' {FormAstromifs},
  uGlobals in 'source\code\uGlobals.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TFormAstromifs, FormAstromifs);
  Application.Run;
end.
