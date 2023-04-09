//
(* PaintedNotes for colorized keyboards https://github.com/geoblock *)
//
program Astromifon;

uses
  Forms,
  uGlobals in 'Source\Code\uGlobals.pas',
  fInitials in 'Source\Interface\fInitials.pas' {FormInitial},
  fAbout in 'source\interface\fAbout.pas' {FormAbouts},
  fDialogs in 'Source\Interface\fDialogs.pas' {FormDialog},
  fOptions in 'Source\Interface\fOptions.pas' {FormOptions},
  fAstrofon in 'source\interface\fAstrofon.pas' {FormScene};

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TFormScene, FormScene);
  Application.Run;
end.
