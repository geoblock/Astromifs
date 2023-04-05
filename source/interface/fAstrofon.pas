//
(* PaintedNotes for colorized keyboards https://github.com/geoblock *)
//
unit fAstrofon;

interface

uses
  Winapi.Windows,
  Winapi.Messages,
  Winapi.MMSystem,
  System.SysUtils,
  System.Variants,
  System.Classes,
  System.ImageList,
  System.Math,


  Vcl.Graphics,
  Vcl.Controls,
  Vcl.Forms,
  Vcl.Dialogs,
  Vcl.ImgList,
  Vcl.ExtCtrls,
  Vcl.ComCtrls,
  Vcl.Menus,
  Vcl.ToolWin,
  Vcl.StdCtrls,
  Vcl.Imaging.jpeg,

  GLS.DCE,
  GLS.SoundManager,
  Sounds.BASS,
  GLS.Collision,
  GLS.FPSMovement,
  Physics.NGDManager,
  GLS.BaseClasses,
  GLS.Cadencer,
  GLS.ParticleFX,
  GLS.Material,

  fAbout,
  fOptions,

  GLS.VectorTypes,
  GLS.VectorGeometry,
  GLS.Scene,
  GLS.SceneViewer,
  GLS.FireFX,
  GLS.Objects,
  GLS.HUDObjects,
  GLS.Graph,
  GLS.VectorFileObjects,
  GLS.GeomObjects,
  GLS.Mesh,
  GLS.Coordinates,
  GLS.SkyDome,
  GLS.Navigator,
  GLS.Keyboard,
  GLS.BitmapFont,
  GLS.WindowsFont,
  GLS.SimpleNavigation,
  GLS.Color,
  GLS.Utils;

type
  TPianoKeySet = set of 0 .. 87;
  TGuitarKeySet = set of 0 .. 149;

type
  TFormScene = class(TForm)
    MainMenu: TMainMenu;
    File1: TMenuItem;
    miNew: TMenuItem;
    miOpen: TMenuItem;
    miSave: TMenuItem;
    miSaveAs: TMenuItem;
    N1: TMenuItem;
    miExit: TMenuItem;
    Edit1: TMenuItem;
    Undo1: TMenuItem;
    N5: TMenuItem;
    Cut1: TMenuItem;
    Copy1: TMenuItem;
    Paste1: TMenuItem;
    N4: TMenuItem;
    GoTo1: TMenuItem;
    N3: TMenuItem;
    Window1: TMenuItem;
    miChangeFPS: TMenuItem;
    miMaxFPS: TMenuItem;
    miMiddleFPS: TMenuItem;
    miMinFPS: TMenuItem;
    N6: TMenuItem;
    Hide1: TMenuItem;
    Show1: TMenuItem;
    Help1: TMenuItem;
    miContents: TMenuItem;
    N2: TMenuItem;
    miAbout: TMenuItem;
    miTools: TMenuItem;
    miCalculate: TMenuItem;
    miOptions: TMenuItem;
    GLSceneViewer1: TGLSceneViewer;
    GLScene: TGLScene;
    LightSource1: TGLLightSource;
    Camera: TGLCamera;
    dcScene: TGLDummyCube;
    spEarth: TGLSphere;
    dcMoon: TGLDummyCube;
    spMoon: TGLSphere;
    GLMatLibChromatic: TGLMaterialLibrary;
    dcPianoKeys: TGLDummyCube;
    cbPianoStand: TGLCube;
    dcCamera: TGLDummyCube;
    GLCadencer: TGLCadencer;
    ControlBar1: TControlBar;
    ToolBarFile: TToolBar;
    PanelBottom: TPanel;
    StatusBar1: TStatusBar;
    PanelLeft: TPanel;
    PanelRight: TPanel;
    Timer: TTimer;
    GLWindowsBitmapFont: TGLWindowsBitmapFont;
    ImageListInterface: TImageList;
    ImageListComponents: TImageList;
    ImageListObjects: TImageList;
    spCore: TGLSphere;
    LightSource2: TGLLightSource;
    TreeView1: TTreeView;
    TreeView2: TTreeView;
    ToolButton1: TToolButton;
    ToolButton2: TToolButton;
    ToolButton3: TToolButton;
    ToolButton4: TToolButton;
    ToolButton5: TToolButton;
    dcGuitarKeys: TGLDummyCube;
    cbGuitarStand: TGLCube;
    grdGuitarKeys: TGLXYZGrid;
    grdPianoKeysH: TGLXYZGrid;
    EdgeXa: TGLArrowLine;
    EdgeYf: TGLArrowLine;
    EdgeZc: TGLArrowLine;
    dcEdges: TGLDummyCube;
    EdgeZcs: TGLArrowLine;
    EdgeXas: TGLArrowLine;
    GLArrowLine5: TGLArrowLine;
    EdgeZd: TGLArrowLine;
    EdgeZds: TGLArrowLine;
    EdgeYfs: TGLArrowLine;
    EdgeYg: TGLArrowLine;
    EdgeYgs: TGLArrowLine;
    EdgeXb: TGLArrowLine;
    EdgeXe: TGLArrowLine;
    GLBitmapFont: TGLBitmapFont;
    GLSkyDome1: TGLSkyDome;
    HUDTextGuitar: TGLHUDText;
    HUDTextPiano: TGLHUDText;
    GLFlatText1: TGLFlatText;
    rgKeyboardColors: TRadioGroup;
    GLMatLibBW: TGLMaterialLibrary;
    grdPianoKeysV: TGLXYZGrid;
    hfBottom: TGLHeightField;
    procedure miAboutClick(Sender: TObject);
    procedure miOptionsClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure GLCadencerProgress(Sender: TObject; const DeltaTime, NewTime: Double);
    procedure GLSceneViewer1MouseDown(Sender: TObject; Button: TMouseButton; Shift: TShiftState;
      X, Y: Integer);
    procedure GLSceneViewer1MouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
    procedure FormMouseWheel(Sender: TObject; Shift: TShiftState; WheelDelta: Integer;
      MousePos: TPoint; var Handled: Boolean);
    procedure TimerTimer(Sender: TObject);
    procedure miExitClick(Sender: TObject);
    procedure GLSceneViewer1MouseUp(Sender: TObject; Button: TMouseButton; Shift: TShiftState; X,
      Y: Integer);
    procedure rgKeyboardColorsClick(Sender: TObject);
  public
    mx: Integer;
    my: Integer;
    MaterialIndex: Integer;
    procedure MakePianoKeys(Sender: TObject);
    procedure MakeGuitarKeys(Sender: TObject);
    property Action;
  private
    CameraHeight: Single;
    PathToAsset: TFileName;

    // piano keyboard
    // 88 keys for 9 octavas or 88 constellations
    NPianoKeys: Integer;
    PianoKeys: array [0 .. 87] of TGLCube;
    KeyDepthW, KeyHeightW, KeyWidthW: Single;
    KeyDepthB, KeyHeightB, KeyWidthB: Single;
    CountOctavas: Integer;

    // guitar fretboard
    // 149 keys for 6 strings
    NGuitarKeys: Integer;
    GuitarKeys: array [0 .. 149] of TGLCube;
    KeyDepth, KeyHeight, KeyWidth: Single;
    StringIndex, FretIndex: Integer;
    CountStrings, CountFrets: Integer;
    OriginX, OriginY, OriginZ: Single;
    OpenMaterials: array [0 .. 5] of Integer; // materials for open strings

    procedure SetPianoKeySizes;
    procedure SetGuitarKeySizes;
  end;

var
  FormScene: TFormScene;
  BlackKeySet: TPianoKeySet;
  NutKeySet: TGuitarKeySet;
  PickDown: TGLCustomSceneObject;
  PickUp: TGLCustomSceneObject;
  old_color, new_color: TVector4f;

const
  SelectionColor: TGLColorVector = (X: 0.243; Y: 0.243; Z: 0.243; W: 1.000);

implementation

{$R *.dfm}

// ----------------------------------------
// Define key sizes ans set for black keys
// ----------------------------------------
procedure TFormScene.SetPianoKeySizes;
begin
  NPianoKeys := 88;

  KeyDepthW := 12;
  KeyDepthB := 6;
  KeyHeightW := 1;
  KeyHeightB := 1;
  KeyWidthW := 1;
  KeyWidthB := 0.6;

  // 36 black keys for 9 octavas
  BlackKeySet :=                                 [{  0 }  1, {  2,   // 0 octava
     3 }  4, {  5 }  6, {  7,  8 }  9, { 10 } 11, { 12 } 13, { 14,   // 1 octava
    15 } 16, { 17 } 18, { 19, 20 } 21, { 22 } 23, { 24 } 25, { 26,   // 2 octava
    27 } 28, { 29 } 30, { 31, 32 } 33, { 34 } 35, { 36 } 37, { 38,   // 3 octava
    39 } 40, { 41 } 42, { 43, 44 } 45, { 46 } 47, { 48 } 49, { 50,   // 4 octava
    51 } 52, { 53 } 54, { 55, 56 } 57, { 58 } 59, { 60 } 61, { 62,   // 5 octava
    63 } 64, { 65 } 66, { 67, 68 } 69, { 70 } 71, { 72 } 73, { 74,   // 6 octava
    75 } 76, { 77 } 78, { 79, 80 } 81, { 82 } 83, { 84 } 85  { 86,   // 7 octava
    87 }];                                                           // 8 octava
end;

// ----------------------------------------------------------------------------------
// Make 88 keys for 3 notes in 0 octave + 84 notes in 7 octavas + 1 note for octava 8
// ----------------------------------------------------------------------------------
procedure TFormScene.MakePianoKeys(Sender: TObject);
var
  i: Integer;
  CurrentX: Single;
begin
  SetPianoKeySizes;

  CurrentX := (cbPianoStand.CubeWidth / 2 - 0.5);
  for i := 0 to NPianoKeys - 1 do
  begin
    PianoKeys[i] := TGLCube.CreateAsChild(dcPianoKeys);
    if (i in BlackKeySet) then
    begin
      PianoKeys[i].Position.SetPoint(CurrentX + 0.5, 1, 3);
      PianoKeys[i].CubeDepth := KeyDepthB;
      PianoKeys[i].CubeHeight := KeyHeightB;
      PianoKeys[i].CubeWidth := KeyWidthB;
      if (rgKeyboardColors.ItemIndex = 1) then
        PianoKeys[i].Material.FrontProperties.Diffuse.Color := clrGray20;
    end
    else
    begin
      PianoKeys[i].Position.SetPoint(CurrentX, 0, 0);
      PianoKeys[i].CubeDepth := KeyDepthW;
      PianoKeys[i].CubeHeight := KeyHeightW;
      PianoKeys[i].CubeWidth := KeyWidthW;
      if (rgKeyboardColors.ItemIndex = 1) then
        PianoKeys[i].Material.FrontProperties.Diffuse.Color := clrWhite;
      CurrentX := CurrentX - KeyWidthW; // moving to the next key position on X axes
    end;
    PianoKeys[i].Tag := i + 1;  // for 1..88 keys
    MaterialIndex := i mod 12;  // or MaterialIndex := (i+3) mod 12;
    case rgKeyboardColors.ItemIndex of
    0: begin
         PianoKeys[i].Material.MaterialLibrary := GLMatLibChromatic;
         PianoKeys[i].Material.LibMaterialName := GLMatLibChromatic.Materials[MaterialIndex].Name;
       end;
    (*
    3: begin
         PianoKeys[i].Material.MaterialLibrary := GLMatLibScriabin;
         PianoKeys[i].Material.LibMaterialName := GLMatLibScriabin.Materials[MaterialIndex].Name;
       end;
    *)
    end;
  end;
end;

// --------------------------------------
// Define guitar key sizes
// --------------------------------------
procedure TFormScene.SetGuitarKeySizes;
begin
  NGuitarKeys := 150;

  KeyDepth := 1;
  KeyHeight := 1;
  KeyWidth := 2;
  OriginX := (cbGuitarStand.CubeWidth / 2 - 1);
  OriginY := 11;
  OriginZ := 5;

  CountStrings := 6;
  CountFrets := 25;

  FretIndex := 0;
  StringIndex := 0;
  MaterialIndex := 0;

  // open strings colors
  OpenMaterials[0] := 7;  // e -> mi
  OpenMaterials[1] := 2;  // b -> si
  OpenMaterials[2] := 10; // g -> sol
  OpenMaterials[3] := 5;  // d -> re
  OpenMaterials[4] := 0;  // a -> la
  OpenMaterials[5] := 7;  // e -> mi
end;

// ---------------------------------------------------
//          Make 150 guitar keys
// ---------------------------------------------------
procedure TFormScene.MakeGuitarKeys(Sender: TObject);
var
  i, j, k, NumString: Integer; // current string number

begin
  SetGuitarKeySizes;
  i := 0;
  j := 0;
  k := 0; // current keyindex
  NumString := 0; // the first string number
  for i := 0 to CountFrets - 1 do
  begin
    for j := 0 to CountStrings - 1 do
    begin
      GuitarKeys[k] := TGLCube.CreateAsChild(dcGuitarKeys);
      NumString := k div 25 + 1;
      FretIndex := k mod 25;
      MaterialIndex := (OpenMaterials[NumString - 1] + FretIndex) mod 12;
      GuitarKeys[k].Material.MaterialLibrary := GLMatLibChromatic;
      GuitarKeys[k].Material.LibMaterialName := GLMatLibChromatic.Materials[MaterialIndex].Name;
      //
      GuitarKeys[k].Position.SetPoint(OriginX - FretIndex * 2 - 1, OriginY - NumString + 1,
        OriginZ);
      GuitarKeys[k].CubeDepth := KeyDepth;
      GuitarKeys[k].CubeHeight := KeyHeight;
      GuitarKeys[k].CubeWidth := KeyWidth;
      Inc(k);
    end;
  end;
  NGuitarKeys := k;
end;

// --------------------------------------------------------------------

procedure TFormScene.FormCreate(Sender: TObject);
begin
  PathToAsset := GetCurrentDir();
  // Loading maps for planets
  SetCurrentDir(PathToAsset + '\asset\map');
  spEarth.Material.Texture.Disabled := False;
  spEarth.Material.Texture.Image.LoadFromFile('earth.jpg');
  spMoon.Material.Texture.Disabled := False;
  spMoon.Material.Texture.Image.LoadFromFile('moon.jpg');

  SetCurrentDir(PathToAsset + '\asset\texture');
  cbPianoStand.Material.Texture.Disabled := False;
  cbPianoStand.Material.Texture.Image.LoadFromFile('ashwood.jpg');
  cbGuitarStand.Material.Texture.Disabled := False;
  cbGuitarStand.Material.Texture.Image.LoadFromFile('ashwood.jpg');

  SetCurrentDir(PathToAsset  + '\asset\font');
  GLBitmapFont.Glyphs.LoadFromFile('darkgold_font.bmp');

  //GLWindowsBitmapFont.Glyphs.LoadFromFile('toonfont.bmp');

  MakePianoKeys(Self);
  MakeGuitarKeys(Self);
  rgKeyboardColorsClick(nil);
end;

// -----------------------------------------------------------------------------------------
procedure TFormScene.FormMouseWheel(Sender: TObject; Shift: TShiftState; WheelDelta: Integer;
  MousePos: TPoint; var Handled: Boolean);
begin
  Camera.AdjustDistanceToTarget(Power(1.1, WheelDelta / 120));
end;

// ----------------------------------------------------------------------------------------
procedure TFormScene.GLCadencerProgress(Sender: TObject; const DeltaTime, NewTime: Double);
var
  speed: Single;
  MyString: String;
begin
  inherited;
  // dcEdges.TurnAngle := -NewTime * 60;
  // spEarth.TurnAngle := newTime * 60;
  if IsKeyDown(VK_SHIFT) then
    speed := 5 * DeltaTime
  else
    speed := DeltaTime;
  with Camera.Position do
  begin // WASD
    if IsKeyDown(68) then
      dcCamera.Translate(Z * speed, 0, -X * speed);
    if IsKeyDown(65) then
      dcCamera.Translate(-Z * speed, 0, X * speed);
    if IsKeyDown(87) then
      dcCamera.Translate(-X * speed, 0, -Z * speed);
    if IsKeyDown(83) then
      dcCamera.Translate(X * speed, 0, Z * speed);
    if IsKeyDown('Q') or IsKeyDown('é') then
    begin
      CameraHeight := CameraHeight + 10 * speed;
      dcCamera.Position.Y := CameraHeight;
    end;
    if IsKeyDown('E') or IsKeyDown('ó') then
    begin
      CameraHeight := CameraHeight - 10 * speed;
      dcCamera.Position.Y := CameraHeight;
    end;
    if IsKeyDown(VK_ESCAPE) then
      Close;
  end;
  StatusBar1.Panels[1].Text := 'X: ' + FloatToStrF(Camera.Position.X, ffFixed, 7, 2);
  StatusBar1.Panels[2].Text := 'Y: ' + FloatToStrF(Camera.Position.Y, ffFixed, 7, 2);
  StatusBar1.Panels[3].Text := 'Z: ' + FloatToStrF(Camera.Position.Z, ffFixed, 7, 2);
end;

// ------------------------------------------------------------------------------

procedure TFormScene.GLSceneViewer1MouseDown(Sender: TObject; Button: TMouseButton;
  Shift: TShiftState; X, Y: Integer);
var
  PickedObject: TGLCustomSceneObject;
begin
  if Button = TMouseButton.mbRight then
  begin
    // if an object is picked...
    PickedObject := (GLSceneViewer1.Buffer.GetPickedObject(X, Y) as TGLCustomSceneObject);
    if Assigned(PickedObject) and (PickedObject is TGLCube) then
    begin
      // ...turn it to yellow and show its name

      PickedObject.Material.FrontProperties.Emission.Color := clrYellow; // not from MatLib

      old_color := PickedObject.Material.FrontProperties.Emission.Color;
      new_color := vectorlerp(old_color, clrWhite, 0.25);

(*
      if Pickdown.Name = 'GLSkyDome1' then
        Pickdown.Material.FrontProperties.Emission.Color := old_color
      else
        Pickdown.Material.FrontProperties.Emission.Color := new_color;
*)
      ShowMessage('You clicked the ' + PickedObject.Name + ' ' + IntToStr(PickedObject.Tag)
        + ' key');
      end;
  end;
  mx := X;
  my := Y;
end;

// ------------------------------------------------------------------------------

procedure TFormScene.GLSceneViewer1MouseUp(Sender: TObject; Button: TMouseButton;
  Shift: TShiftState; X, Y: Integer);
begin
  if Button = TMouseButton.mbRight then // check if the mouse button is still pressed
  begin
    PickUp := (GLSceneViewer1.Buffer.GetPickedObject(X, Y) as TGLCustomSceneObject);
    if Assigned(PickUp) then
      if (PickUp = PickDown) and (PickUp is TGLCube) then
      begin
        SetCurrentDir(PathToAsset + '\assets\audio\piano');

        PickUp.Material.FrontProperties.Emission.Color := old_color;
        if PickUp.Name = 'cbPianoStand' then
        begin
//          S1 := S1 + '1';
          // sndPlaySound('red.wav',snd_ASync);
          PlaySound('k30.ogg', hinstance, SND_RESOURCE or SND_ASYNC);
        end;
      end;
  end;
end;

// ------------------------------------------------------------------------------

procedure TFormScene.GLSceneViewer1MouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
begin
  if (ssLeft in Shift) then
    Camera.MoveAroundTarget(my - Y, mx - X);
  mx := X;
  my := Y;
end;

// ------------------------------------------------------------------------------

procedure TFormScene.miExitClick(Sender: TObject);
begin
  FormScene.Close;
end;

// ------------------------------------------------------------------------------

procedure TFormScene.miOptionsClick(Sender: TObject);
begin
  inherited;
  with TFormOptions.Create(nil) do
    try
      ShowModal;
    finally
      Free;
    end;
end;

// ------------------------------------------------------------------------------

procedure TFormScene.miAboutClick(Sender: TObject);
begin
  with TFormAbouts.Create(nil) do
    try
      ShowModal;
    finally
      Free;
    end;
end;

// ------------------------------------------------------------------------------

procedure TFormScene.rgKeyboardColorsClick(Sender: TObject);
begin
  if rgKeyboardColors.ItemIndex <> 2 then
  begin
    grdPianoKeysH.LineColor.Color := clrBlack;
    grdPianoKeysV.LineColor.Color := clrBlack
  end
  else
  begin
    grdPianoKeysH.LineColor.Color := clrWhite;;
    grdPianoKeysV.LineColor.Color := clrWhite;
  end;

  MakePianoKeys(Sender);
  GLSceneViewer1.Invalidate;
end;


// ------------------------------------------------------------------------------

procedure TFormScene.TimerTimer(Sender: TObject);
begin
  StatusBar1.Panels[0].Text := Format('FPS:  %.1f ', [GLSceneViewer1.FramesPerSecond]);
  // Format('%d particles, %.1f FPS', [GLParticles1.Count, GLSceneViewer1.FramesPerSecond]);
  GLSceneViewer1.ResetPerformanceMonitor;
end;

end.
