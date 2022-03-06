//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 25 14:20:47 2017 by ROOT version 6.04/14
// from TTree output/output
// found on file: /xrootdfs/atlas/dq2/rucio/user/hwang43/f8/23/user.hwang43.11724839._000001.hist-output.root
//////////////////////////////////////////////////////////

#ifndef CalibrationAnalysis_h
#define CalibrationAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class CalibrationAnalysis {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Fixed size dimensions of array or collections stored in the TTree if any.

        // Declaration of leaf types
        Int_t           m_mcChannelNumber;
        Int_t           runNumber;
        ULong64_t       eventNumber;
        Double_t        m_mgg;
        Double_t        m_met;
        Double_t        m_weight;
        Double_t        m_sumet;
        Double_t        truth_ph_pt1;
        Double_t        truth_ph_eta1;
        Double_t        truth_ph_phi1;
        Double_t        truth_ph_e1;
        Double_t        truth_ph_pt2;
        Double_t        truth_ph_eta2;
        Double_t        truth_ph_phi2;
        Double_t        truth_ph_e2;
        Int_t           n_ph;
        Double_t        ph_t1;
        Double_t        ph_t2;
        Double_t        ph_z1;
        Double_t        ph_z2;
        Double_t        ph_pt1;
        Double_t        ph_pt2;
        Double_t        ph_eta1;
        Double_t        ph_eta2;
        Double_t        ph_phi1;
        Double_t        ph_phi2;
        Int_t           ph_convType1;
        Int_t           ph_convType2;
        Bool_t          ph_isTight1;
        Bool_t          ph_isTight2;
        Bool_t          ph_isLoose1;
        Bool_t          ph_isLoose2;
        Bool_t          ph_passIso1;
        Bool_t          ph_passIso2;
        Double_t        ph_conv_z1;
        Double_t        ph_conv_z_err1;
        Double_t        ph_calo_z1;
        Double_t        ph_calo_z_err1;
        Double_t        ph_conv_z2;
        Double_t        ph_conv_z_err2;
        Double_t        ph_calo_z2;
        Double_t        ph_calo_z_err2;
        Float_t         ph1_topoetcone20;
        Float_t         ph1_topoetcone40;
        Float_t         ph1_ptcone20;
        Float_t         ph2_topoetcone20;
        Float_t         ph2_topoetcone40;
        Float_t         ph2_ptcone20;
        Float_t         ph1_f1;
        Float_t         ph1_f3;
        Float_t         ph1_etas2;
        Float_t         ph1_phis2;
        ULong64_t       ph1_maxEcell_onlId;
        Int_t           ph1_maxEcell_gain;
        Float_t         ph1_maxEcell_E;
        Float_t         ph1_maxEcell_x;
        Float_t         ph1_maxEcell_y;
        Float_t         ph1_maxEcell_z;
        Float_t         ph1_maxEcell_time;
        Float_t         ph2_f1;
        Float_t         ph2_f3;
        Float_t         ph2_etas2;
        Float_t         ph2_phis2;
        ULong64_t       ph2_maxEcell_onlId;
        Int_t           ph2_maxEcell_gain;
        Float_t         ph2_maxEcell_E;
        Float_t         ph2_maxEcell_x;
        Float_t         ph2_maxEcell_y;
        Float_t         ph2_maxEcell_z;
        Float_t         ph2_maxEcell_time;
        Double_t        PV_z;
        Double_t        NPV;
        Double_t        wt_xs;
        Double_t        wt_kf;
        Double_t        wt_ge;
        Double_t        wt_mc;
        Double_t        wt_nEvents;
        Double_t        wt_wt;
        Double_t        intlumi;
        Int_t           HLT_2g20_tight;
        Int_t           HLT_2g50_loose;
        Int_t           HLT_g100_loose;
        Int_t           HLT_g10_loose;
        Int_t           HLT_g120_loose;
        Int_t           HLT_g140_loose;
        Int_t           HLT_g15_loose_L1EM7;
        Int_t           HLT_g20_loose_L1EM12;
        Int_t           HLT_g25_loose_L1EM15;
        Int_t           HLT_g35_loose_L1EM15;
        Int_t           HLT_g35_loose_g25_loose;
        Int_t           HLT_g35_medium_g25_medium;
        Int_t           HLT_g40_loose_L1EM15;
        Int_t           HLT_g45_loose_L1EM15;
        Int_t           HLT_g50_loose_L1EM15;
        Int_t           HLT_g60_loose;
        Int_t           HLT_g70_loose;
        Int_t           HLT_g80_loose;

        // List of branches
        TBranch        *b_m_mcChannelNumber;   //!
        TBranch        *b_runNumber;   //!
        TBranch        *b_eventNumber;   //!
        TBranch        *b_m_mgg;   //!
        TBranch        *b_m_met;   //!
        TBranch        *b_m_weight;   //!
        TBranch        *b_m_sumet;   //!
        TBranch        *b_truth_ph_pt1;   //!
        TBranch        *b_truth_ph_eta1;   //!
        TBranch        *b_truth_ph_phi1;   //!
        TBranch        *b_truth_ph_e1;   //!
        TBranch        *b_truth_ph_pt2;   //!
        TBranch        *b_truth_ph_eta2;   //!
        TBranch        *b_truth_ph_phi2;   //!
        TBranch        *b_truth_ph_e2;   //!
        TBranch        *b_n_ph;   //!
        TBranch        *b_ph_t1;   //!
        TBranch        *b_ph_t2;   //!
        TBranch        *b_ph_z1;   //!
        TBranch        *b_ph_z2;   //!
        TBranch        *b_ph_pt1;   //!
        TBranch        *b_ph_pt2;   //!
        TBranch        *b_ph_eta1;   //!
        TBranch        *b_ph_eta2;   //!
        TBranch        *b_ph_phi1;   //!
        TBranch        *b_ph_phi2;   //!
        TBranch        *b_ph_convType1;   //!
        TBranch        *b_ph_convType2;   //!
        TBranch        *b_ph_isTight1;   //!
        TBranch        *b_ph_isTight2;   //!
        TBranch        *b_ph_isLoose1;   //!
        TBranch        *b_ph_isLoose2;   //!
        TBranch        *b_ph_passIso1;   //!
        TBranch        *b_ph_passIso2;   //!
        TBranch        *b_ph_conv_z1;   //!
        TBranch        *b_ph_conv_z_err1;   //!
        TBranch        *b_ph_calo_z1;   //!
        TBranch        *b_ph_calo_z_err1;   //!
        TBranch        *b_ph_conv_z2;   //!
        TBranch        *b_ph_conv_z_err2;   //!
        TBranch        *b_ph_calo_z2;   //!
        TBranch        *b_ph_calo_z_err2;   //!
        TBranch        *b_ph1_topoetcone20;   //!
        TBranch        *b_ph1_topoetcone40;   //!
        TBranch        *b_ph1_ptcone20;   //!
        TBranch        *b_ph2_topoetcone20;   //!
        TBranch        *b_ph2_topoetcone40;   //!
        TBranch        *b_ph2_ptcone20;   //!
        TBranch        *b_ph1_f1;   //!
        TBranch        *b_ph1_f3;   //!
        TBranch        *b_ph1_etas2;   //!
        TBranch        *b_ph1_phis2;   //!
        TBranch        *b_ph1_maxEcell_onlId;   //!
        TBranch        *b_ph1_maxEcell_gain;   //!
        TBranch        *b_ph1_maxEcell_E;   //!
        TBranch        *b_ph1_maxEcell_x;   //!
        TBranch        *b_ph1_maxEcell_y;   //!
        TBranch        *b_ph1_maxEcell_z;   //!
        TBranch        *b_ph1_maxEcell_time;   //!
        TBranch        *b_ph2_f1;   //!
        TBranch        *b_ph2_f3;   //!
        TBranch        *b_ph2_etas2;   //!
        TBranch        *b_ph2_phis2;   //!
        TBranch        *b_ph2_maxEcell_onlId;   //!
        TBranch        *b_ph2_maxEcell_gain;   //!
        TBranch        *b_ph2_maxEcell_E;   //!
        TBranch        *b_ph2_maxEcell_x;   //!
        TBranch        *b_ph2_maxEcell_y;   //!
        TBranch        *b_ph2_maxEcell_z;   //!
        TBranch        *b_ph2_maxEcell_time;   //!
        TBranch        *b_PV_z;   //!
        TBranch        *b_NPV;   //!
        TBranch        *b_wt_xs;   //!
        TBranch        *b_wt_kf;   //!
        TBranch        *b_wt_ge;   //!
        TBranch        *b_wt_mc;   //!
        TBranch        *b_wt_nEvents;   //!
        TBranch        *b_wt_wt;   //!
        TBranch        *b_intlumi;   //!
        TBranch        *b_HLT_2g20_tight;   //!
        TBranch        *b_HLT_2g50_loose;   //!
        TBranch        *b_HLT_g100_loose;   //!
        TBranch        *b_HLT_g10_loose;   //!
        TBranch        *b_HLT_g120_loose;   //!
        TBranch        *b_HLT_g140_loose;   //!
        TBranch        *b_HLT_g15_loose_L1EM7;   //!
        TBranch        *b_HLT_g20_loose_L1EM12;   //!
        TBranch        *b_HLT_g25_loose_L1EM15;   //!
        TBranch        *b_HLT_g35_loose_L1EM15;   //!
        TBranch        *b_HLT_g35_loose_g25_loose;   //!
        TBranch        *b_HLT_g35_medium_g25_medium;   //!
        TBranch        *b_HLT_g40_loose_L1EM15;   //!
        TBranch        *b_HLT_g45_loose_L1EM15;   //!
        TBranch        *b_HLT_g50_loose_L1EM15;   //!
        TBranch        *b_HLT_g60_loose;   //!
        TBranch        *b_HLT_g70_loose;   //!
        TBranch        *b_HLT_g80_loose;   //!

        CalibrationAnalysis(TTree *tree=0);
        virtual ~CalibrationAnalysis();
        virtual Int_t    Cut(Long64_t entry);
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(TTree *tree);
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show(Long64_t entry = -1);




        //HISTOGRAM INITIALIZATION

        TH1F* h_t_uncorr; //!
        TH1F* h_t_pass0; //!
        TH1F* h_t_pass0_fromFile; //!
        TH1F* h_t_pass1; //!
        TH2F* h_t_pass1_vs_e; //!
        TH2F* h_t_pass1_vs_et; //!
        TProfile* p_t_pass0_vs_runN; //!
        TProfile* p_t_pass1_vs_runN; //!
        TH1F* h_t_pass0_emb; //!
        TH1F* h_t_pass1_emb; //!
        TH1F* h_t_pass0_emec; //!
        TH1F* h_t_pass1_emec; //!
        TProfile* p_t_pass0_vs_runN_emb; //!
        TProfile* p_t_pass1_vs_runN_emb; //!
        TProfile* p_t_pass0_vs_runN_emec; //!
        TProfile* p_t_pass1_vs_runN_emec; //!
        TH1F* h_t_pass0_emba; //!
        TH1F* h_t_pass1_emba; //!
        TH1F* h_t_pass0_emeca; //!
        TH1F* h_t_pass1_emeca; //!
        TProfile* p_t_pass0_vs_runN_emba; //!
        TProfile* p_t_pass1_vs_runN_emba; //!
        TProfile* p_t_pass0_vs_runN_emeca; //!
        TProfile* p_t_pass1_vs_runN_emeca; //!
        TH1F* h_t_pass0_embc; //!
        TH1F* h_t_pass1_embc; //!
        TH1F* h_t_pass0_emecc; //!
        TH1F* h_t_pass1_emecc; //!
        TProfile* p_t_pass0_vs_runN_embc; //!
        TProfile* p_t_pass1_vs_runN_embc; //!
        TProfile* p_t_pass0_vs_runN_emecc; //!
        TProfile* p_t_pass1_vs_runN_emecc; //!


        TProfile* p_t_pass1_vs_e_emba; //!
        TProfile* p_t_pass1_vs_e_embc; //!
        TProfile* p_t_pass1_vs_e_emeca; //!
        TProfile* p_t_pass1_vs_e_emecc; //!

        TProfile* p_t_pass1_vs_febN; //!
        TProfile* p_t_pass2_vs_febN; //!
        TH1F* h_t_pass2; //!


        TH1F* h_t_pass1_val; //!
        TProfile* p_t_pass1_vs_febN_val; //!
        TProfile* p_t_pass2_vs_febN_val; //!
        TH1F* h_t_pass2_val; //!


        TH1F* h_t_pass2_narrow; //!
        TH1F* h_t_pass2_wide; //!


        TH1F* h_t_uncorr_val; //!


        TProfile* p_t_pass1_vs_runN_val; //!

        TH1F* h_t_pass1_emb_val; //!

        TH1F* h_t_pass1_emec_val; //!

        TProfile* p_t_pass1_vs_runN_emb_val; //!

        TProfile* p_t_pass1_vs_runN_emec_val; //!

        TH1F* h_t_pass1_emba_val; //!

        TH1F* h_t_pass1_emeca_val; //!

        TProfile* p_t_pass1_vs_runN_emba_val; //!

        TProfile* p_t_pass1_vs_runN_emeca_val; //!

        TH1F* h_t_pass1_embc_val; //!

        TH1F* h_t_pass1_emecc_val; //!

        TProfile* p_t_pass1_vs_runN_embc_val; //!

        TProfile* p_t_pass1_vs_runN_emecc_val; //!


        TH1F* h_t; //!

        //Good timing
        TH1F* h_t_good; //!
        TH1F* h_e_good; //!
        TH1F* h_et_good; //!
        TH1F* h_met_good; //!
        TH1F* h_Ecell_E_good; //!
        TH1F* h_Ecell_x_good; //!
        TH1F* h_Ecell_y_good; //!
        TH1F* h_Ecell_z_good; //!
        TH1F* h_eta_good; //!
        TH1F* h_phi_good; //!
        TH1F* h_cl_eta_good; //!
        TH1F* h_cl_phi_good; //!
        TH1F* h_f1_good; //!
        TH1F* h_f3_good; //!


        //Bad timing
        TH1F* h_t_bad; //!
        TH1F* h_e_bad; //!
        TH1F* h_et_bad; //!
        TH1F* h_met_bad; //!
        TH1F* h_Ecell_E_bad; //!
        TH1F* h_Ecell_x_bad; //!
        TH1F* h_Ecell_y_bad; //!
        TH1F* h_Ecell_z_bad; //!
        TH1F* h_eta_bad; //!
        TH1F* h_phi_bad; //!
        TH1F* h_cl_eta_bad; //!
        TH1F* h_cl_phi_bad; //!
        TH1F* h_f1_bad; //!
        TH1F* h_f3_bad; //!



        TH1F* h_t_bad_stdm2; //!




        TH1F* h_t_32febs; //!
        TH1F* h_t_corr_32febs; //!


        TH1F* h_t_32febs_wide; //!
        TH1F* h_t_corr_32febs_wide; //!


        TH1F* h_t_fn[256]; //!

        TH1F* h_t_febs_slot[8]; //!
        TH1F* h_t_corr_febs_slot[8]; //!
        TH1F* h_e_slot[8]; //!
        TH1F* h_e_cell_slot[8]; //!
        TH1F* h_et_slot[8]; //!
        TH1F* h_eta_slot[8]; //!
        TH1F* h_phi_slot[8]; //!
        TProfile* p_t_vs_E_slot[8]; //!
        TProfile* p_t_corr_vs_E_slot[8]; //!
        TProfile* p_t_vs_ET_slot[8]; //!
        TProfile* p_t_corr_vs_ET_slot[8]; //!
        TProfile* p_t_vs_E_cell_slot[8]; //!
        TProfile* p_t_pass2_vs_E_cell_slot[8]; //!
        TH1F* h_t_pass3_slot[8]; //!
        TProfile* p_t_pass3_vs_E_cell_slot[8]; //!
        TProfile* p_t_pass3_vs_channel_no_slot[8]; //!
        TH1F* h_t_pass4_slot[8]; //!
        TProfile* p_t_pass4_vs_channel_no_slot[8]; //!
        TH1F* h_t_pass3_slot_val[8]; //!
        TProfile* p_t_pass3_vs_E_cell_slot_val[8]; //!
        TProfile* p_t_pass3_vs_channel_no_slot_val[8]; //!
        TH1F* h_t_pass4_slot_val[8]; //!
        TProfile* p_t_pass4_vs_channel_no_slot_val[8]; //!

        TProfile* p_t_pass3_vs_channel_no_sym_slot[4]; //!
        TProfile* p_t_pass4_vs_channel_no_sym_slot[4]; //!
        TProfile* p_t_pass3_vs_channel_no_sym_slot_val[4]; //!
        TProfile* p_t_pass4_vs_channel_no_sym_slot_val[4]; //!



        TH1F* h_t_pass3; //!
        TH1F* h_t_pass3_emb; //!
        TH1F* h_t_pass4; //!
        TH1F* h_t_pass4_emb; //!

        TH1F* h_t_pass3_val; //!
        TH1F* h_t_pass3_emb_val; //!
        TH1F* h_t_pass4_val; //!
        TH1F* h_t_pass4_emb_val; //!

        TProfile* p_t_pass3_vs_del_eta; //!
        TProfile* p_t_pass3_vs_del_phi; //!
        TProfile* p_t_pass3_vs_f1; //!
        TProfile* p_t_pass3_vs_f3; //!

        TProfile* p_t_pass4_vs_del_eta; //!
        TProfile* p_t_pass4_vs_del_phi; //!
        TProfile* p_t_pass4_vs_f1; //!
        TProfile* p_t_pass4_vs_f3; //!





        TH1F* h_t_febs_narrow; //!
        TH1F* h_t_corr_febs_narrow; //!
        TH1F* h_e_narrow; //!
        TH1F* h_e_cell_narrow; //!
        TH1F* h_et_narrow; //!
        TH1F* h_eta_narrow; //!
        TH1F* h_phi_narrow; //!

        TH1F* h_t_febs_wide; //!
        TH1F* h_t_corr_febs_wide; //!
        TH1F* h_e_wide; //!
        TH1F* h_e_cell_wide; //!
        TH1F* h_et_wide; //!
        TH1F* h_eta_wide; //!
        TH1F* h_phi_wide; //!

        TH1F* h_t_febs_emec; //!
        TH1F* h_t_corr_febs_emec; //!
        TH1F* h_e_emec; //!
        TH1F* h_e_cell_emec; //!
        TH1F* h_et_emec; //!
        TH1F* h_eta_emec; //!
        TH1F* h_phi_emec; //!
        TProfile* p_t_corr_emec_vs_subslot; //!

        TH1F* h_e_fullFEB; //!






        //USER-DEFINED FUNCTIONS
        void bookHists(); //!
        void writeHists(); //!
        bool badOnlidCheck(unsigned int onl_id); //!
        double tofCorrection(double x, double y, double z, double t_uncorr, double pv_z); //!
        double getLineFromFile(string fileName, int l_no);
        int getFebNo(unsigned long int onlId);
        int getSlotNo(unsigned long int onlId);
        int getSubslotNo(unsigned long int onlId);
        int getChannelNo(unsigned long int onlId);



};

#endif

#ifdef CalibrationAnalysis_cxx
CalibrationAnalysis::CalibrationAnalysis(TTree *tree) : fChain(0) 
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        //      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Ntuples/exot6.root");
        //      if (!f || !f->IsOpen()) {
        TFile *f = new TFile("/data/users/akahn/ntuples/exot6.root");
        //      }
        f->GetObject("output",tree);

    }
    Init(tree);
}

CalibrationAnalysis::~CalibrationAnalysis()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t CalibrationAnalysis::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t CalibrationAnalysis::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void CalibrationAnalysis::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("m_mcChannelNumber", &m_mcChannelNumber, &b_m_mcChannelNumber);
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    fChain->SetBranchAddress("m_mgg", &m_mgg, &b_m_mgg);
    fChain->SetBranchAddress("m_met", &m_met, &b_m_met);
    fChain->SetBranchAddress("m_weight", &m_weight, &b_m_weight);
    fChain->SetBranchAddress("m_sumet", &m_sumet, &b_m_sumet);
    fChain->SetBranchAddress("truth_ph_pt1", &truth_ph_pt1, &b_truth_ph_pt1);
    fChain->SetBranchAddress("truth_ph_eta1", &truth_ph_eta1, &b_truth_ph_eta1);
    fChain->SetBranchAddress("truth_ph_phi1", &truth_ph_phi1, &b_truth_ph_phi1);
    fChain->SetBranchAddress("truth_ph_e1", &truth_ph_e1, &b_truth_ph_e1);
    fChain->SetBranchAddress("truth_ph_pt2", &truth_ph_pt2, &b_truth_ph_pt2);
    fChain->SetBranchAddress("truth_ph_eta2", &truth_ph_eta2, &b_truth_ph_eta2);
    fChain->SetBranchAddress("truth_ph_phi2", &truth_ph_phi2, &b_truth_ph_phi2);
    fChain->SetBranchAddress("truth_ph_e2", &truth_ph_e2, &b_truth_ph_e2);
    fChain->SetBranchAddress("n_ph", &n_ph, &b_n_ph);
    fChain->SetBranchAddress("ph_t1", &ph_t1, &b_ph_t1);
    fChain->SetBranchAddress("ph_t2", &ph_t2, &b_ph_t2);
    fChain->SetBranchAddress("ph_z1", &ph_z1, &b_ph_z1);
    fChain->SetBranchAddress("ph_z2", &ph_z2, &b_ph_z2);
    fChain->SetBranchAddress("ph_pt1", &ph_pt1, &b_ph_pt1);
    fChain->SetBranchAddress("ph_pt2", &ph_pt2, &b_ph_pt2);
    fChain->SetBranchAddress("ph_eta1", &ph_eta1, &b_ph_eta1);
    fChain->SetBranchAddress("ph_eta2", &ph_eta2, &b_ph_eta2);
    fChain->SetBranchAddress("ph_phi1", &ph_phi1, &b_ph_phi1);
    fChain->SetBranchAddress("ph_phi2", &ph_phi2, &b_ph_phi2);
    fChain->SetBranchAddress("ph_convType1", &ph_convType1, &b_ph_convType1);
    fChain->SetBranchAddress("ph_convType2", &ph_convType2, &b_ph_convType2);
    fChain->SetBranchAddress("ph_isTight1", &ph_isTight1, &b_ph_isTight1);
    fChain->SetBranchAddress("ph_isTight2", &ph_isTight2, &b_ph_isTight2);
    fChain->SetBranchAddress("ph_isLoose1", &ph_isLoose1, &b_ph_isLoose1);
    fChain->SetBranchAddress("ph_isLoose2", &ph_isLoose2, &b_ph_isLoose2);
    fChain->SetBranchAddress("ph_passIso1", &ph_passIso1, &b_ph_passIso1);
    fChain->SetBranchAddress("ph_passIso2", &ph_passIso2, &b_ph_passIso2);
    fChain->SetBranchAddress("ph_conv_z1", &ph_conv_z1, &b_ph_conv_z1);
    fChain->SetBranchAddress("ph_conv_z_err1", &ph_conv_z_err1, &b_ph_conv_z_err1);
    fChain->SetBranchAddress("ph_calo_z1", &ph_calo_z1, &b_ph_calo_z1);
    fChain->SetBranchAddress("ph_calo_z_err1", &ph_calo_z_err1, &b_ph_calo_z_err1);
    fChain->SetBranchAddress("ph_conv_z2", &ph_conv_z2, &b_ph_conv_z2);
    fChain->SetBranchAddress("ph_conv_z_err2", &ph_conv_z_err2, &b_ph_conv_z_err2);
    fChain->SetBranchAddress("ph_calo_z2", &ph_calo_z2, &b_ph_calo_z2);
    fChain->SetBranchAddress("ph_calo_z_err2", &ph_calo_z_err2, &b_ph_calo_z_err2);
    fChain->SetBranchAddress("ph1_topoetcone20", &ph1_topoetcone20, &b_ph1_topoetcone20);
    fChain->SetBranchAddress("ph1_topoetcone40", &ph1_topoetcone40, &b_ph1_topoetcone40);
    fChain->SetBranchAddress("ph1_ptcone20", &ph1_ptcone20, &b_ph1_ptcone20);
    fChain->SetBranchAddress("ph2_topoetcone20", &ph2_topoetcone20, &b_ph2_topoetcone20);
    fChain->SetBranchAddress("ph2_topoetcone40", &ph2_topoetcone40, &b_ph2_topoetcone40);
    fChain->SetBranchAddress("ph2_ptcone20", &ph2_ptcone20, &b_ph2_ptcone20);
    fChain->SetBranchAddress("ph1_f1", &ph1_f1, &b_ph1_f1);
    fChain->SetBranchAddress("ph1_f3", &ph1_f3, &b_ph1_f3);
    fChain->SetBranchAddress("ph1_etas2", &ph1_etas2, &b_ph1_etas2);
    fChain->SetBranchAddress("ph1_phis2", &ph1_phis2, &b_ph1_phis2);
    fChain->SetBranchAddress("ph1_maxEcell_onlId", &ph1_maxEcell_onlId, &b_ph1_maxEcell_onlId);
    fChain->SetBranchAddress("ph1_maxEcell_gain", &ph1_maxEcell_gain, &b_ph1_maxEcell_gain);
    fChain->SetBranchAddress("ph1_maxEcell_E", &ph1_maxEcell_E, &b_ph1_maxEcell_E);
    fChain->SetBranchAddress("ph1_maxEcell_x", &ph1_maxEcell_x, &b_ph1_maxEcell_x);
    fChain->SetBranchAddress("ph1_maxEcell_y", &ph1_maxEcell_y, &b_ph1_maxEcell_y);
    fChain->SetBranchAddress("ph1_maxEcell_z", &ph1_maxEcell_z, &b_ph1_maxEcell_z);
    fChain->SetBranchAddress("ph1_maxEcell_time", &ph1_maxEcell_time, &b_ph1_maxEcell_time);
    fChain->SetBranchAddress("ph2_f1", &ph2_f1, &b_ph2_f1);
    fChain->SetBranchAddress("ph2_f3", &ph2_f3, &b_ph2_f3);
    fChain->SetBranchAddress("ph2_etas2", &ph2_etas2, &b_ph2_etas2);
    fChain->SetBranchAddress("ph2_phis2", &ph2_phis2, &b_ph2_phis2);
    fChain->SetBranchAddress("ph2_maxEcell_onlId", &ph2_maxEcell_onlId, &b_ph2_maxEcell_onlId);
    fChain->SetBranchAddress("ph2_maxEcell_gain", &ph2_maxEcell_gain, &b_ph2_maxEcell_gain);
    fChain->SetBranchAddress("ph2_maxEcell_E", &ph2_maxEcell_E, &b_ph2_maxEcell_E);
    fChain->SetBranchAddress("ph2_maxEcell_x", &ph2_maxEcell_x, &b_ph2_maxEcell_x);
    fChain->SetBranchAddress("ph2_maxEcell_y", &ph2_maxEcell_y, &b_ph2_maxEcell_y);
    fChain->SetBranchAddress("ph2_maxEcell_z", &ph2_maxEcell_z, &b_ph2_maxEcell_z);
    fChain->SetBranchAddress("ph2_maxEcell_time", &ph2_maxEcell_time, &b_ph2_maxEcell_time);
    fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
    fChain->SetBranchAddress("NPV", &NPV, &b_NPV);
    fChain->SetBranchAddress("wt_xs", &wt_xs, &b_wt_xs);
    fChain->SetBranchAddress("wt_kf", &wt_kf, &b_wt_kf);
    fChain->SetBranchAddress("wt_ge", &wt_ge, &b_wt_ge);
    fChain->SetBranchAddress("wt_mc", &wt_mc, &b_wt_mc);
    fChain->SetBranchAddress("wt_nEvents", &wt_nEvents, &b_wt_nEvents);
    fChain->SetBranchAddress("wt_wt", &wt_wt, &b_wt_wt);
    fChain->SetBranchAddress("intlumi", &intlumi, &b_intlumi);
    fChain->SetBranchAddress("HLT_2g20_tight", &HLT_2g20_tight, &b_HLT_2g20_tight);
    fChain->SetBranchAddress("HLT_2g50_loose", &HLT_2g50_loose, &b_HLT_2g50_loose);
    fChain->SetBranchAddress("HLT_g100_loose", &HLT_g100_loose, &b_HLT_g100_loose);
    fChain->SetBranchAddress("HLT_g10_loose", &HLT_g10_loose, &b_HLT_g10_loose);
    fChain->SetBranchAddress("HLT_g120_loose", &HLT_g120_loose, &b_HLT_g120_loose);
    fChain->SetBranchAddress("HLT_g140_loose", &HLT_g140_loose, &b_HLT_g140_loose);
    fChain->SetBranchAddress("HLT_g15_loose_L1EM7", &HLT_g15_loose_L1EM7, &b_HLT_g15_loose_L1EM7);
    fChain->SetBranchAddress("HLT_g20_loose_L1EM12", &HLT_g20_loose_L1EM12, &b_HLT_g20_loose_L1EM12);
    fChain->SetBranchAddress("HLT_g25_loose_L1EM15", &HLT_g25_loose_L1EM15, &b_HLT_g25_loose_L1EM15);
    fChain->SetBranchAddress("HLT_g35_loose_L1EM15", &HLT_g35_loose_L1EM15, &b_HLT_g35_loose_L1EM15);
    fChain->SetBranchAddress("HLT_g35_loose_g25_loose", &HLT_g35_loose_g25_loose, &b_HLT_g35_loose_g25_loose);
    fChain->SetBranchAddress("HLT_g35_medium_g25_medium", &HLT_g35_medium_g25_medium, &b_HLT_g35_medium_g25_medium);
    fChain->SetBranchAddress("HLT_g40_loose_L1EM15", &HLT_g40_loose_L1EM15, &b_HLT_g40_loose_L1EM15);
    fChain->SetBranchAddress("HLT_g45_loose_L1EM15", &HLT_g45_loose_L1EM15, &b_HLT_g45_loose_L1EM15);
    fChain->SetBranchAddress("HLT_g50_loose_L1EM15", &HLT_g50_loose_L1EM15, &b_HLT_g50_loose_L1EM15);
    fChain->SetBranchAddress("HLT_g60_loose", &HLT_g60_loose, &b_HLT_g60_loose);
    fChain->SetBranchAddress("HLT_g70_loose", &HLT_g70_loose, &b_HLT_g70_loose);
    fChain->SetBranchAddress("HLT_g80_loose", &HLT_g80_loose, &b_HLT_g80_loose);
    Notify();
}

Bool_t CalibrationAnalysis::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void CalibrationAnalysis::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t CalibrationAnalysis::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef CalibrationAnalysis_cxx
