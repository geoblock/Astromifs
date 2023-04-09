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
  fSonofon;

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
    GLMatLib: TGLMaterialLibrary;
    PanelTop: TPanel;
    PanelBottom: TPanel;
    RadioGroup1: TRadioGroup;
    grbShow: TGroupBox;
    chbStars: TCheckBox;
    CheckBox1: TCheckBox;
    CheckBox2: TCheckBox;
    sbMilkiWay: TGLSkyBox;
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
    miAstrofon: TMenuItem;
    Hide2: TMenuItem;
    Show2: TMenuItem;
    N7: TMenuItem;
    procedure miAboutClick(Sender: TObject);
    procedure Open1Click(Sender: TObject);
    procedure Save1Click(Sender: TObject);
    procedure SaveAs1Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure GLCadencerProgress(Sender: TObject; const DeltaTime, NewTime: Double);
    procedure tvConstellationsClick(Sender: TObject);
    procedure miAstrofonClick(Sender: TObject);
  private
    CurrentPath: TFileName;
    PathToData: TFileName;
    function LoadTexture(Matname, Filename: string): TGLLibMaterial;
  public
    procedure HandleKeys(d: Double);
  end;

var
  frmAstroD: TfrmAstroD;

implementation

{$R *.dfm}

//-----------------------------------------------------------------------

function TfrmAstroD.LoadTexture(Matname, Filename: string): TGLLibMaterial;
begin
  Result := GLMatLib.AddTextureMaterial(Matname, Filename);
  Result.Material.Texture.Disabled := False;
  Result.Material.Texture.TextureMode := tmDecal;
end;

//-----------------------------------------------------------------------

procedure TfrmAstroD.FormCreate(Sender: TObject);
var
  SkyColture, PlanetMap: TFileName;
begin
  PathToData := GetCurrentDir() + '\data';
  SetCurrentDir(PathToData);

  // Skybox cubemaps
  GLMatLib.TexturePaths := GetCurrentDir() + '\cubemap';
  LoadTexture('Left', 'mw_left.jpg');
  LoadTexture('Right', 'mw_right.jpg');
  LoadTexture('Top', 'mw_top.jpg');
  LoadTexture('Bottom', 'mw_bot.jpg');
  LoadTexture('Front', 'mw_front.jpg');
  LoadTexture('Back', 'mw_back.jpg');

  PlanetMap := PathToData + '\map\mars.jpg';

  Planet.Material.Texture.Disabled := False;
  Planet.Material.Texture.Image.LoadFromFile(PlanetMap);

  SkyColture := PathToData + '\constellation\ConstellationNames.dat';
  tvConstellations.LoadFromFile(SkyColture);

  SkyColture := PathToData + '\constellation\ConstShortNames.dat';
  tvCurrent.LoadFromFile(SkyColture);

end;

//-----------------------------------------------------------------------

procedure TfrmAstroD.Open1Click(Sender: TObject);
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

procedure TfrmAstroD.tvConstellationsClick(Sender: TObject);
begin
  //
end;

//-----------------------------------------------------------------------

procedure TfrmAstroD.Save1Click(Sender: TObject);
begin
  // Save TreeView
end;

//-----------------------------------------------------------------------

procedure TfrmAstroD.SaveAs1Click(Sender: TObject);
begin
  // Save TreeView As...
end;

//-----------------------------------------------------------------------


//-----------------------------------------------------------------------

procedure TfrmAstroD.GLCadencerProgress(Sender: TObject; const DeltaTime, NewTime: Double);
begin
 //
  HandleKeys(deltaTime);
  GLUserInterface1.Mouselook;
  GLUserInterface1.MouseUpdate;
  GLSceneViewer.Invalidate;

end;

//-----------------------------------------------------------------------

procedure TfrmAstroD.HandleKeys(d: Double);
begin
  if IsKeyDown('W') or IsKeyDown('Z') then
    Camera.Move(d);
  if IsKeyDown('S') then
    Camera.Move(-d);
  if IsKeyDown('A') or IsKeyDown('A') then
    Camera.Slide(-d);
  if IsKeyDown('D') then
    Camera.Slide(d);

  if IsKeyDown(VK_ESCAPE) then
    Close;
end;


//----------------------------------------------------------------

procedure TfrmAstroD.miAstrofonClick(Sender: TObject);
begin
  with TfrmAstrofon.Create(Self) do
    try
      ShowModal;
    finally
      Free;
    end;
end;

//-----------------------------------------------------------------------

procedure TfrmAstroD.miAboutClick(Sender: TObject);
begin
  // Revived constellations from myths
  with TfrmAbout.Create(nil) do
    try
      ShowModal;
    finally
      Free;
    end;
end;

end.
