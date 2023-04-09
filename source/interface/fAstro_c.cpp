//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "fAstro_c.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "GLS.BaseClasses"
#pragma link "GLS.Cadencer"
#pragma link "GLS.Material"
#pragma link "GLS.Scene"
#pragma link "GLS.SceneViewer"
#pragma resource "*.dfm"
TfrmAstroC *frmAstroC;
//---------------------------------------------------------------------------
__fastcall TfrmAstroC::TfrmAstroC(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
