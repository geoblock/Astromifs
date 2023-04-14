object FormAstromifs: TFormAstromifs
  Left = 0
  Top = 0
  Caption = 'Astromyths & Legends'
  ClientHeight = 622
  ClientWidth = 906
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -12
  Font.Name = 'Segoe UI'
  Font.Style = []
  Menu = MainMenu1
  Position = poScreenCenter
  OnCreate = FormCreate
  TextHeight = 15
  object PanelLeft: TPanel
    Left = 0
    Top = 44
    Width = 161
    Height = 559
    Align = alLeft
    TabOrder = 0
    object tvConstellations: TTreeView
      Left = 1
      Top = 42
      Width = 159
      Height = 475
      Align = alClient
      Indent = 19
      TabOrder = 0
      OnClick = tvConstellationsClick
    end
    object PanelTop: TPanel
      Left = 1
      Top = 1
      Width = 159
      Height = 41
      Align = alTop
      Caption = 'Constellations'
      TabOrder = 1
    end
    object PanelBottom: TPanel
      Left = 1
      Top = 517
      Width = 159
      Height = 41
      Align = alBottom
      TabOrder = 2
    end
  end
  object ControlBar1: TControlBar
    Left = 0
    Top = 0
    Width = 906
    Height = 44
    Align = alTop
    TabOrder = 1
  end
  object StatusBar1: TStatusBar
    Left = 0
    Top = 603
    Width = 906
    Height = 19
    Panels = <>
  end
  object PanelRight: TPanel
    Left = 752
    Top = 44
    Width = 154
    Height = 559
    Align = alRight
    TabOrder = 3
    object RadioGroup1: TRadioGroup
      Left = 16
      Top = 272
      Width = 113
      Height = 127
      Caption = 'Coordinate grid'
      ItemIndex = 0
      Items.Strings = (
        'Equatorial'
        'Ecliptic'
        'Galactic'
        'Supergalactic')
      TabOrder = 0
    end
    object grbShow: TGroupBox
      Left = 16
      Top = 405
      Width = 113
      Height = 121
      Caption = 'View'
      TabOrder = 1
      object chbStars: TCheckBox
        Left = 13
        Top = 24
        Width = 97
        Height = 17
        Caption = 'Stars'
        TabOrder = 0
      end
      object CheckBox1: TCheckBox
        Left = 13
        Top = 56
        Width = 97
        Height = 17
        Caption = 'Lines'
        TabOrder = 1
      end
      object CheckBox2: TCheckBox
        Left = 13
        Top = 88
        Width = 97
        Height = 17
        Caption = 'Borders'
        TabOrder = 2
      end
    end
    object tvCurrent: TTreeView
      Left = 1
      Top = 42
      Width = 152
      Height = 224
      Align = alTop
      Indent = 19
      TabOrder = 2
    end
    object PanelTopR: TPanel
      Left = 1
      Top = 1
      Width = 152
      Height = 41
      Align = alTop
      Caption = 'Current'
      TabOrder = 3
    end
  end
  object GLSceneViewer: TGLSceneViewer
    Left = 161
    Top = 44
    Width = 591
    Height = 559
    Camera = Camera
    Buffer.BackgroundColor = clBlack
    FieldOfView = 155.768493652343800000
    PenAsTouch = False
    Align = alClient
    TabOrder = 4
  end
  object MainMenu1: TMainMenu
    Left = 232
    Top = 40
    object File1: TMenuItem
      Caption = '&File'
      object New1: TMenuItem
        Caption = '&New'
      end
      object Open1: TMenuItem
        Caption = '&Open...'
        OnClick = Open1Click
      end
      object Save1: TMenuItem
        Caption = '&Save'
        OnClick = Save1Click
      end
      object SaveAs1: TMenuItem
        Caption = 'Save &As...'
        OnClick = SaveAs1Click
      end
      object N2: TMenuItem
        Caption = '-'
      end
      object Print1: TMenuItem
        Caption = '&Print...'
      end
      object PrintSetup1: TMenuItem
        Caption = 'P&rint Setup...'
      end
      object N1: TMenuItem
        Caption = '-'
      end
      object Exit1: TMenuItem
        Caption = 'E&xit'
        OnClick = Exit1Click
      end
    end
    object Edit1: TMenuItem
      Caption = '&Edit'
      object Undo1: TMenuItem
        Caption = '&Undo'
        ShortCut = 16474
      end
      object Repeat1: TMenuItem
        Caption = '&Repeat <command>'
      end
      object N5: TMenuItem
        Caption = '-'
      end
      object Cut1: TMenuItem
        Caption = 'Cu&t'
        ShortCut = 16472
      end
      object Copy1: TMenuItem
        Caption = '&Copy'
        ShortCut = 16451
      end
      object Paste1: TMenuItem
        Caption = '&Paste'
        ShortCut = 16470
      end
      object PasteSpecial1: TMenuItem
        Caption = 'Paste &Special...'
      end
      object N4: TMenuItem
        Caption = '-'
      end
      object Find1: TMenuItem
        Caption = '&Find...'
      end
      object Replace1: TMenuItem
        Caption = 'R&eplace...'
      end
      object GoTo1: TMenuItem
        Caption = '&Go To...'
      end
      object N3: TMenuItem
        Caption = '-'
      end
      object Links1: TMenuItem
        Caption = 'Lin&ks...'
      end
      object Object1: TMenuItem
        Caption = '&Object'
      end
    end
    object Window2: TMenuItem
      Caption = '&View'
      object miSonofon: TMenuItem
        Caption = 'Sonofon...'
        OnClick = miSonofonClick
      end
      object N7: TMenuItem
        Caption = '-'
      end
      object Hide2: TMenuItem
        Caption = '&Hide'
      end
      object Show2: TMenuItem
        Caption = '&Show...'
      end
    end
    object Window1: TMenuItem
      Caption = '&Window'
      object NewWindow1: TMenuItem
        Caption = '&New Window'
      end
      object Tile1: TMenuItem
        Caption = '&Tile'
      end
      object Cascade1: TMenuItem
        Caption = '&Cascade'
      end
      object ArrangeAll1: TMenuItem
        Caption = '&Arrange All'
      end
      object N6: TMenuItem
        Caption = '-'
      end
      object Hide1: TMenuItem
        Caption = '&Hide'
      end
      object Show1: TMenuItem
        Caption = '&Show...'
      end
    end
    object Help2: TMenuItem
      Caption = '&Help'
      object Contents2: TMenuItem
        Caption = '&Contents'
      end
      object Tutorial1: TMenuItem
        Caption = '&Tutorial'
      end
      object miAbout: TMenuItem
        Caption = '&About...'
        OnClick = miAboutClick
      end
    end
  end
  object GLScene: TGLScene
    Left = 178
    Top = 64
    object SkyDome: TGLSkyDome
      Visible = False
      Bands = <
        item
          StartColor.Color = {0000803F0000803F0000803F0000803F}
          StopAngle = 15.000000000000000000
        end
        item
          StartAngle = 15.000000000000000000
          StopAngle = 90.000000000000000000
          StopColor.Color = {938C0C3E938C0C3E938E0E3F0000803F}
          Stacks = 4
        end>
      Stars = <>
    end
    object SkyBoxMilkyWay: TGLSkyBox
      MaterialLibrary = GLMatLib
      MatNameTop = 'Top'
      MatNameBottom = 'Bottom'
      MatNameLeft = 'Left'
      MatNameRight = 'Right'
      MatNameFront = 'Front'
      MatNameBack = 'Back'
      MatNameClouds = 'Clouds'
      CloudsPlaneOffset = 0.200000002980232200
      CloudsPlaneSize = 32.000000000000000000
    end
    object Camera: TGLCamera
      DepthOfView = 100000.000000000000000000
      FocalLength = 60.000000000000000000
      NearPlaneBias = 0.100000001490116100
      TargetObject = dcWorld
      Position.Coordinates = {0000004000000000000000000000803F}
    end
    object LightSource: TGLLightSource
      ConstAttenuation = 1.000000000000000000
      Position.Coordinates = {0000C8420000C8420000C8420000803F}
      SpotCutOff = 180.000000000000000000
      object LensFlare: TGLLensFlare
        Seed = 1465
        FlareIsNotOccluded = True
      end
    end
    object dcWorld: TGLDummyCube
      CubeSize = 1.000000000000000000
      object Planet: TGLSphere
        Radius = 0.500000000000000000
        Slices = 64
        Stacks = 64
      end
      object ffPlanet: TGLFreeForm
        Visible = False
      end
    end
  end
  object GLCadencer: TGLCadencer
    Scene = GLScene
    OnProgress = GLCadencerProgress
    Left = 296
    Top = 200
  end
  object GLMatLib: TGLMaterialLibrary
    Left = 352
    Top = 340
  end
  object GLNavigator: TGLNavigator
    Left = 448
    Top = 56
  end
  object GLUserInterface1: TGLUserInterface
    GLNavigator = GLNavigator
    Left = 584
    Top = 56
  end
  object GLSimpleNavigation1: TGLSimpleNavigation
    Form = Owner
    GLSceneViewer = GLSceneViewer
    FormCaption = 'Astromyths & Legends - %FPS'
    KeyCombinations = <
      item
        ShiftState = [ssLeft, ssRight]
        Action = snaZoom
      end
      item
        ShiftState = [ssLeft]
        Action = snaMoveAroundTarget
      end
      item
        ShiftState = [ssRight]
        Action = snaMoveAroundTarget
      end>
    Left = 330
    Top = 426
  end
  object OpenDialog: TOpenDialog
    Left = 195
    Top = 312
  end
end
