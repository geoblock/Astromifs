//---------------------------------------------------------------------------

#ifndef fAstro_xH
#define fAstro_xH
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <FMX.Controls.hpp>
#include <FMX.Forms.hpp>
#include <FMX.Controls.Presentation.hpp>
#include <FMX.Layouts.hpp>
#include <FMX.Menus.hpp>
#include <FMX.StdCtrls.hpp>
#include <FMX.TreeView.hpp>
#include <FMX.Types.hpp>
//---------------------------------------------------------------------------
class TfrmAstroX : public TForm
{
__published:	// IDE-managed Components
	TMainMenu *MainMenu1;
	TPanel *PanelLeft;
	TPanel *PanelRight;
	TTreeView *TreeView1;
private:	// User declarations
public:		// User declarations
	__fastcall TfrmAstroX(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TfrmAstroX *frmAstroX;
//---------------------------------------------------------------------------
#endif
