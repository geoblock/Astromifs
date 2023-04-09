//---------------------------------------------------------------------------

#ifndef fAstro_cH
#define fAstro_cH
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include "GLS.BaseClasses.hpp"
#include "GLS.Cadencer.hpp"
#include "GLS.Material.hpp"
#include "GLS.Scene.hpp"
#include "GLS.SceneViewer.hpp"
#include <Vcl.ComCtrls.hpp>
#include <Vcl.ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TfrmAstroC : public TForm
{
__published:	// IDE-managed Components
	TGLSceneViewer *GLSceneViewer1;
	TGLScene *GLScene1;
	TGLMaterialLibrary *GLMaterialLibrary1;
	TPanel *PanelLeft;
	TPanel *PanelRight;
	TStatusBar *StatusBar1;
	TControlBar *ControlBar1;
	TGLCadencer *GLCadencer1;
private:	// User declarations
public:		// User declarations
	__fastcall TfrmAstroC(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmAstroC *frmAstroC;
//---------------------------------------------------------------------------
#endif
