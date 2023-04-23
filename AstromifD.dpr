program AstromifD;

uses
  Vcl.Forms,
  fAstromifs in 'source\interface\fAstromifs.pas' {FormAstromifs},
  uGlobals in 'source\code\uGlobals.pas',
  uConstellations in 'source\code\uConstellations.pas',
  fAbout in 'source\interface\fAbout.pas' {frmAbout};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TFormAstromifs, FormAstromifs);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
