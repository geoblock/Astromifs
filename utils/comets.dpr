program comets;

uses
  Vcl.Forms,
  comets_f in 'comets_f.pas' {frmComets};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmComets, frmComets);
  Application.Run;
end.
