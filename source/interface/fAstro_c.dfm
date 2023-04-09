object frmAstroC: TfrmAstroC
  Left = 0
  Top = 0
  Caption = 'Astromyths & Legends C'
  ClientHeight = 602
  ClientWidth = 830
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -12
  Font.Name = 'Segoe UI'
  Font.Style = []
  TextHeight = 15
  object GLSceneViewer1: TGLSceneViewer
    Left = 137
    Top = 33
    Width = 551
    Height = 550
    Buffer.BackgroundColor = clBlack
    PenAsTouch = False
    Align = alClient
    TabOrder = 0
  end
  object PanelLeft: TPanel
    Left = 0
    Top = 33
    Width = 137
    Height = 550
    Align = alLeft
    TabOrder = 1
    ExplicitHeight = 381
  end
  object PanelRight: TPanel
    Left = 688
    Top = 33
    Width = 142
    Height = 550
    Align = alRight
    TabOrder = 2
    ExplicitLeft = 480
    ExplicitHeight = 381
  end
  object StatusBar1: TStatusBar
    Left = 0
    Top = 583
    Width = 830
    Height = 19
    Panels = <>
    ExplicitLeft = 232
    ExplicitTop = 384
    ExplicitWidth = 0
  end
  object ControlBar1: TControlBar
    Left = 0
    Top = 0
    Width = 830
    Height = 33
    Align = alTop
    TabOrder = 4
    ExplicitWidth = 622
  end
  object GLScene1: TGLScene
    Left = 176
    Top = 120
  end
  object GLMaterialLibrary1: TGLMaterialLibrary
    Left = 200
    Top = 216
  end
  object GLCadencer1: TGLCadencer
    Left = 336
    Top = 136
  end
end
