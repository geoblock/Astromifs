program Astros;

uses
  Vcl.Forms,
  fAstros in 'source\vcl\fAstros.pas' {frmAstros},
  Astro.Matlib in 'source\Astro.Matlib.pas',
  Astro.Sphlib in 'source\Astro.Sphlib.pas',
  Astro.Pnulib in 'source\Astro.Pnulib.pas',
  Astro.Delib in 'source\Astro.Delib.pas',
  Astro.Plalib in 'source\Astro.Plalib.pas',
  Astro.Keplib in 'source\Astro.Keplib.pas',
  Astro.Moonlib in 'source\Astro.Moonlib.pas',
  Astro.Sunlib in 'source\Astro.Sunlib.pas',
  Astro.P15lib in 'source\Astro.P15lib.pas',
  Astro.P69lib in 'source\Astro.P69lib.pas',
  Astro.Phylib in 'source\Astro.Phylib.pas',
  Astro.Timlib in 'source\Astro.Timlib.pas',
  fAbout in 'source\vcl\fAbout.pas' {frmAbout};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmAstros, frmAstros);
  Application.CreateForm(TfrmAbout, frmAbout);
  Application.Run;
end.
