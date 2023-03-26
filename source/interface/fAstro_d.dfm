object frmAstroD: TfrmAstroD
  Left = 0
  Top = 0
  Caption = 'Astromyths & Legends D'
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
      Items.NodeData = {
        0301000000300000000000000000000000FFFFFFFFFFFFFFFF00000000000000
        0000000000010941006E00640072006F006D00650064006100}
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
      Top = 26
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
      Top = 192
      Width = 113
      Height = 169
      Caption = 'View'
      TabOrder = 1
      object chbStars: TCheckBox
        Left = 16
        Top = 24
        Width = 97
        Height = 17
        Caption = 'Stars'
        TabOrder = 0
      end
      object CheckBox1: TCheckBox
        Left = 16
        Top = 56
        Width = 97
        Height = 17
        Caption = 'Lines'
        TabOrder = 1
      end
      object CheckBox2: TCheckBox
        Left = 16
        Top = 88
        Width = 97
        Height = 17
        Caption = 'Borders'
        TabOrder = 2
      end
    end
  end
  object GLSceneViewer1: TGLSceneViewer
    Left = 161
    Top = 44
    Width = 591
    Height = 559
    Buffer.BackgroundColor = clBlack
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
      object Index1: TMenuItem
        Caption = '&Index'
      end
      object Commands1: TMenuItem
        Caption = 'Co&mmands'
      end
      object Procedures1: TMenuItem
        Caption = '&Procedures'
      end
      object Keyboard1: TMenuItem
        Caption = '&Keyboard'
      end
      object Tutorial1: TMenuItem
        Caption = '&Tutorial'
      end
      object HowtoUseHelp2: TMenuItem
        Caption = '&How to Use Help'
      end
      object About2: TMenuItem
        Caption = '&About...'
        OnClick = About2Click
      end
    end
  end
  object GLScene1: TGLScene
    Left = 296
    Top = 128
  end
  object GLCadencer1: TGLCadencer
    Left = 296
    Top = 200
  end
  object GLMaterialLibrary1: TGLMaterialLibrary
    Left = 288
    Top = 272
  end
end
