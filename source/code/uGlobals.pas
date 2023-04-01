//
(* PaintedNotes for colorized keyboards https://github.com/geoblock *)
//
unit uGlobals;

interface

uses
  System.SysUtils,
  System.IniFiles;


const
  //cRegistryKey = 'Software\SoundCube';
  RegSoundCube = PathDelim + 'SOFTWARE' + PathDelim + 'Gramofon' + PathDelim;

  SoundCubeSize = 100;

var
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
  ftSQL = 'sql';  // sql
  ftTVN = 'tvn';  // treeview nodes

//==========================================================================
implementation
//==========================================================================



//---------------------------
 initialization

   FormatSettings.DecimalSeparator := '.';


end.
