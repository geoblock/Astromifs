unit fAstro_d;

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

  Astro.Matlib,
  Vcl.Menus,
  Vcl.ComCtrls,
  Vcl.ExtCtrls,
  GLS.Material,
  GLS.Cadencer,
  GLS.BaseClasses,
  GLS.Scene,
  GLS.SceneViewer,
  Vcl.StdCtrls;

type
  TfrmAstroD = class(TForm)
    MainMenu1: TMainMenu;
    File1: TMenuItem;
    New1: TMenuItem;
    Open1: TMenuItem;
    Save1: TMenuItem;
    SaveAs1: TMenuItem;
    Print1: TMenuItem;
    PrintSetup1: TMenuItem;
    Exit1: TMenuItem;
    N1: TMenuItem;
    N2: TMenuItem;
    Edit1: TMenuItem;
    Undo1: TMenuItem;
    Repeat1: TMenuItem;
    Cut1: TMenuItem;
    Copy1: TMenuItem;
    Paste1: TMenuItem;
    PasteSpecial1: TMenuItem;
    Find1: TMenuItem;
    Replace1: TMenuItem;
    GoTo1: TMenuItem;
    Links1: TMenuItem;
    Object1: TMenuItem;
    N3: TMenuItem;
    N4: TMenuItem;
    N5: TMenuItem;
    Window1: TMenuItem;
    NewWindow1: TMenuItem;
    Tile1: TMenuItem;
    Cascade1: TMenuItem;
    ArrangeAll1: TMenuItem;
    Hide1: TMenuItem;
    Show1: TMenuItem;
    N6: TMenuItem;
    Help2: TMenuItem;
    Contents2: TMenuItem;
    Index1: TMenuItem;
    Commands1: TMenuItem;
    Procedures1: TMenuItem;
    Keyboard1: TMenuItem;
    Tutorial1: TMenuItem;
    HowtoUseHelp2: TMenuItem;
    About2: TMenuItem;
    PanelLeft: TPanel;
    tvConstellations: TTreeView;
    ControlBar1: TControlBar;
    StatusBar1: TStatusBar;
    PanelRight: TPanel;
    GLSceneViewer1: TGLSceneViewer;
    GLScene1: TGLScene;
    GLCadencer1: TGLCadencer;
    GLMaterialLibrary1: TGLMaterialLibrary;
    PanelTop: TPanel;
    PanelBottom: TPanel;
    RadioGroup1: TRadioGroup;
    grbShow: TGroupBox;
    chbStars: TCheckBox;
    CheckBox1: TCheckBox;
    CheckBox2: TCheckBox;
    procedure About2Click(Sender: TObject);
    procedure Open1Click(Sender: TObject);
    procedure Save1Click(Sender: TObject);
    procedure SaveAs1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  frmAstroD: TfrmAstroD;

implementation

{$R *.dfm}

procedure TfrmAstroD.About2Click(Sender: TObject);
begin
  // Revived constellations from myths
end;

procedure TfrmAstroD.Open1Click(Sender: TObject);
begin
  // Load TreeView ...
end;

procedure TfrmAstroD.Save1Click(Sender: TObject);
begin
  // Save TreeView
end;

procedure TfrmAstroD.SaveAs1Click(Sender: TObject);
begin
  // Save TreeView As...
end;

end.
