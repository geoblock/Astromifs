unit fAstromifs;

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

  Vcl.Menus,
  Vcl.ComCtrls,
  Vcl.ExtCtrls,
  GLS.Material,
  GLS.Cadencer,
  GLS.BaseClasses,
  GLS.Scene,
  GLS.SceneViewer,
  Vcl.StdCtrls,

  GLS.Keyboard,
  GLS.Coordinates,
  GLS.Texture,
  GLS.SkyDome,
  GLS.Navigator,
  GLS.LensFlare,
  GLS.Objects,
  GLS.SimpleNavigation,

  fAbout,
  uGlobals,
  GLS.VectorFileObjects;

type
  TFormAstromifs = class(TForm)
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
    Tutorial1: TMenuItem;
    miAbout: TMenuItem;
    PanelLeft: TPanel;
    tvConstellations: TTreeView;
    ControlBar1: TControlBar;
    StatusBar1: TStatusBar;
    PanelRight: TPanel;
    GLSceneViewer: TGLSceneViewer;
    GLScene: TGLScene;
    GLCadencer: TGLCadencer;
    PanelTop: TPanel;
    PanelBottom: TPanel;
    RadioGroup1: TRadioGroup;
    grbShow: TGroupBox;
    chbStars: TCheckBox;
    CheckBox1: TCheckBox;
    CheckBox2: TCheckBox;
    GLNavigator: TGLNavigator;
    GLUserInterface1: TGLUserInterface;
    Camera: TGLCamera;
    LightSource: TGLLightSource;
    dcWorld: TGLDummyCube;
    LensFlare: TGLLensFlare;
    Planet: TGLSphere;
    GLSimpleNavigation1: TGLSimpleNavigation;
    OpenDialog: TOpenDialog;
    tvCurrent: TTreeView;
    PanelTopR: TPanel;
    Window2: TMenuItem;
    Hide2: TMenuItem;
    Show2: TMenuItem;
    N7: TMenuItem;
    SkyDome: TGLSkyDome;
    ffPlanet: TGLFreeForm;
    ConstellationLines: TGLLines;
    ConstellationBorders: TGLLines;
    procedure miAboutClick(Sender: TObject);
    procedure Open1Click(Sender: TObject);
    procedure Save1Click(Sender: TObject);
    procedure SaveAs1Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure GLCadencerProgress(Sender: TObject; const DeltaTime, NewTime: Double);
    procedure tvConstellationsClick(Sender: TObject);
    procedure Exit1Click(Sender: TObject);
  private
  public
    procedure HandleKeys(d: Double);
  end;

var
  FormAstromifs: TFormAstromifs;

implementation

{$R *.dfm}

//-----------------------------------------------------------------------

procedure TFormAstromifs.Exit1Click(Sender: TObject);
begin
  Close;
end;

procedure TFormAstromifs.FormCreate(Sender: TObject);
var
  ConstNames, PlanetMap: TFileName;
begin
  PathToData := GetCurrentDir() + '\data';
  CurrentPath := PathToData;
  SetCurrentDir(CurrentPath + '\cubemap');

  // Skybox stars
  SkyDome.Visible := True;
  SkyDome.Bands.Clear;

  PlanetMap := CurrentPath + '\map\earth.jpg';

  Planet.Material.Texture.Disabled := False;
  Planet.Material.Texture.Image.LoadFromFile(PlanetMap);

  //if FileExists(FileName) then
//    SkyDome.Stars.LoadStarsFile(FileName);
  Catalog := PathToData + '\catalog\hipparcos.stars';
  if FileExists(Catalog) then
  begin
    SkyDome.Bands.Clear;
    SkyDome.Stars.Clear;
    SkyDome.Stars.LoadStarsFile(Catalog);
    SkyDome.StructureChanged;
  end;


  ConstNames := CurrentPath + '\constellation\ConstNames.dat';
    tvConstellations.LoadFromFile(ConstNames);
  ConstNames := CurrentPath + '\constellation\ConstShortNames.dat';
    tvCurrent.LoadFromFile(ConstNames);

  ffPlanet.Assign(Planet);

end;

//-----------------------------------------------------------------------

procedure TFormAstromifs.Open1Click(Sender: TObject);
begin
  // Load next skyculture for constellations ...
  OpenDialog.Filter := 'Constellation (*.dat)|*.dat';
  OpenDialog.InitialDir := PathToData;
  OpenDialog.DefaultExt := '*.dat';
  if OpenDialog.Execute then
  begin
    tvConstellations.LoadFromFile(OpenDialog.FileName);
    CurrentPath := ExtractFilePath(OpenDialog.FileName);
    tvConstellations.Select(tvConstellations.Items[0]);  // goto to new mif
    tvConstellationsClick(Sender);
  end;
end;

//-----------------------------------------------------------------------

procedure TFormAstromifs.tvConstellationsClick(Sender: TObject);
begin
  //
end;

//-----------------------------------------------------------------------

procedure TFormAstromifs.Save1Click(Sender: TObject);
begin
  // Save TreeView
end;

//-----------------------------------------------------------------------

procedure TFormAstromifs.SaveAs1Click(Sender: TObject);
begin
  // Save TreeView As...
end;

//-----------------------------------------------------------------------


//-----------------------------------------------------------------------

procedure TFormAstromifs.GLCadencerProgress(Sender: TObject; const DeltaTime, NewTime: Double);
begin
 //
  HandleKeys(deltaTime);
  GLUserInterface1.Mouselook;
  GLUserInterface1.MouseUpdate;
  GLSceneViewer.Invalidate;

end;

//-----------------------------------------------------------------------

procedure TFormAstromifs.HandleKeys(d: Double);
begin
  if (IsKeyDown('W') or IsKeyDown('Z')) then
    Camera.Move(d);
  if (IsKeyDown('S')) then
    Camera.Move(-d);
  if (IsKeyDown('A') or IsKeyDown('A')) then
    Camera.Slide(-d);
  if (IsKeyDown('D')) then
    Camera.Slide(d);

  if IsKeyDown(VK_ESCAPE) then
    Close;
end;


//-----------------------------------------------------------------------

procedure TFormAstromifs.miAboutClick(Sender: TObject);
begin
  // Revived constellations from myths
  with TfrmAbout.Create(nil) do
    try
      ///LabelTitle.Caption := 'Astromifs';
      ShowModal;
    finally
      Free;
    end;
end;

end.
