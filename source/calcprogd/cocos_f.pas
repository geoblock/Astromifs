unit cocos_f;

interface

uses
  Winapi.Windows,
  Winapi.Messages,
  System.SysUtils,
  System.Variants,
  System.Classes,
  Vcl.Graphics,
  Vcl.Controls,
  Vcl.Forms,
  Vcl.Dialogs,
  Vcl.StdCtrls,

  Astro.Matlib,
  Astro.Pnulib,
  Astro.Sphlib,
  Astro.Sunlib,
  Astro.Timlib;

var
      X,Y,Z,XS,YS,ZS: Double;
      T,TEQX,TEQXN  : Double;
      LS,BS,RS      : Double;
      A             : Double33;
      ECLIPT        : Boolean;
      MODE          : CHAR;


type
  TfrmCocos = class(TForm)
    Edit1: TEdit;
    Memo1: TMemo;
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  frmCocos: TfrmCocos;

implementation

{$R *.dfm}

end.
