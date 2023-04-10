//
(* PaintedNotes for colorized keyboards https://github.com/geoblock *)
//
unit uGlobals;

interface

uses
  System.SysUtils,
  System.IniFiles;


const
  //cRegistryKey = 'Software\Astromifs';
  RegSoundCube = PathDelim + 'SOFTWARE' + PathDelim + 'Astromifs' + PathDelim;

  CubeSize = 1000;

var
  PathToData : TFileName;
  ExePath: TFileName;
  ModelPath:   TFileName;
  TexturePath: TFileName;
  IniFile: TIniFile;

  Language: integer;
  GeneralSection: string = RegSoundCube + 'General';

  SplashStart : Boolean;
  TipOfTheDay : Boolean;

const
  // file types to import/export
  ftTXT = 'txt';  // text
  ftCSV = 'csv';  // csv
  ftDAT = 'dat';  // dat
  ftSQL = 'sqlite';  // sqlite
  ftTVN = 'tvn';  // treeview nodes

//==========================================================================
implementation
//==========================================================================



//---------------------------
 initialization

   FormatSettings.DecimalSeparator := '.';


end.
