program cocos;

uses
  Vcl.Forms,
  cocos_f in 'cocos_f.pas' {frmCocos};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmCocos, frmCocos);
  Application.Run;
end.
