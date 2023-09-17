// Load script that links all the required libraries.
// Adapted from the installation of elSpectro
// [https://github.com/dglazier/elSpectro]
// by Derek Glazier and others
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@alumni.iu.edu
// -----------------------------------------------------------------------------

void Load()
{
    TString lib_ext   = gSystem->GetSoExt();

    //----------------------------------------------------------------------
    // Core library

    TString main_dir  = gSystem->Getenv("ANALYTICRT");

    // Load the main library files
    TString lib  = main_dir + "/lib/libANALYTICRT." + lib_ext;

    // Headers
    TString core = main_dir + "/src/core"; 
    TString phys = main_dir + "/models";
    TString data = main_dir + "/data";

    if (!gSystem->AccessPathName(lib.Data()))
    {
        Int_t lib_loaded = gSystem->Load(lib.Data());
        if (lib_loaded < 0) Fatal("analyticRT::Load", "Library not loaded sucessfully!");

        gInterpreter->AddIncludePath( core.Data());
        gInterpreter->AddIncludePath( phys.Data());
        gInterpreter->AddIncludePath( data.Data());
    }
    else
    {
        Warning("analyticRT::Load", "analyticRT library not found! Looked in: %s", lib.Data());
    }
}