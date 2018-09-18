#ifndef slic3r_MainFrame_hpp_
#define slic3r_MainFrame_hpp_

#include <string>
#include <vector>
#include "PrintConfig.hpp"
#include "../../libslic3r/Utils.hpp"
#include "GUI.hpp"

#include <wx/frame.h>

#include <map>

// class wxApp;
// class wxWindow;
// class wxFrame;
class wxMenuBar;
class wxNotebook;
class wxPanel;
class wxMenu;
class wxProgressDialog;
// class wxComboCtrl;
// class wxString;
// class wxArrayString;
// class wxArrayLong;
// class wxColour;
// class wxBoxSizer;
// class wxFlexGridSizer;
// class wxButton;
// class wxFileDialog;
// class wxStaticBitmap;
// class wxFont;

namespace Slic3r {

class ProgressStatusBar;
class AppController;
//     class AppConfig;
//     class PresetUpdater;
//     class DynamicPrintConfig;
//     class TabIface;

#define _(s)    Slic3r::GUI::I18N::translate((s))

namespace GUI
{
enum QuickSlice
{
    qsUndef,
    qsReslice,
    qsSaveAs,
    qsExportSVG,
    qsExportPNG
};

class MainFrame : public wxFrame
{
    bool        m_no_plater;
    bool        m_loaded;
    int         m_lang_ch_event;
    int         m_preferences_event;

    wxPanel*    m_plater = nullptr;
    wxString    m_qs_last_input_file = wxEmptyString;
    wxString    m_qs_last_output_file = wxEmptyString;
    wxString    m_last_config = wxEmptyString;
    
//     AppConfig*  m_config;

    ProgressStatusBar*              m_statusbar;
    AppController*                  m_appController = nullptr;
    std::map<std::string, Tab*>     m_options_tabs;

    wxMenuItem* append_menu_item(wxMenu* menu,
                                 int id,
                                 const wxString& string,
                                 const wxString& description,
                                 std::function<void(wxCommandEvent& event)> cb,
                                 const std::string& icon = "");

    wxMenuItem* m_menu_item_reslice_now = nullptr;
    wxMenu*     m_plater_menu = nullptr;
    wxMenu*     m_object_menu = nullptr;
    wxMenu*     m_viewMenu = nullptr;

    std::string     get_base_name(const wxString full_name) const ;
    std::string     get_dir_name(const wxString full_name) const ;
public:
    MainFrame() {}
    MainFrame(  bool    no_plater,
                bool    loaded,
                int     lang_ch_event,
                int     preferences_event);
    ~MainFrame() {}


    void        init_tabpanel();
    void        init_menubar();

    void        update_ui_from_settings();
    bool        is_loaded() const { return m_loaded; }
    bool        is_last_input_file() const  { return !m_qs_last_input_file.IsEmpty(); }

    void        on_plater_selection_changed(const bool have_selection);
    void        slice_to_png();
    void        quick_slice(const /*QuickSlice*/int qs = qsUndef/*const bool reslice,
                            const bool export_svg,
                            const bool export_png,
                            const bool save_as*/);
    void        reslice_now();
    void        repair_stl();
    void        export_config();
    void        load_config_file(wxString file = wxEmptyString);
    void        export_configbundle();
    void        load_configbundle(wxString file = wxEmptyString);
    void        load_config(const DynamicPrintConfig& config);
    void        select_tab(size_t tab) const;
    void        select_view(const std::string& direction);


    wxNotebook*         m_tabpanel = nullptr;

    wxProgressDialog*   m_progress_dialog = nullptr;
};

}
}

#endif // slic3r_MainFrame_hpp_