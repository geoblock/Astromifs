program Qosmox;

uses
  System.StartUpCopy,
  FMX.Forms,
  fQosmox in 'source\fmx\fQosmox.pas' {Form1};

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
