#define TimingAnalysis_cxx
#include "Header/TimingAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <string.h>


void TimingAnalysis(string fileName){
    class TimingAnalysis t(0, fileName.c_str());
    t.Loop(fileName);
}


void TimingAnalysis::Loop(string fileName)
{
    //   In a ROOT session, you can do:
    //      root> .L RunTimingNtuple.C
    //      root> RunTimingNtuple t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0) return;


    bool emb = false;
    bool emec = false;
    bool logain = false;
    bool medgain = false;
    bool higain = false;
    bool looseIso = false;
    bool tightIso = false;
    bool tightCaloIso = false;
    bool isData = true;

    int nTot = 0;
    int nLoose = 0;
    int nTightCalo = 0;
    int nTight = 0;

    //MED GAIN CALIBRATION
    TProfile* p_t_vs_et = new TProfile("p_t_vs_e", "Profile of t_{#gamma}^{corrected} vs E, Medium Gain, EMB, Loose Photons", 100, 0., 2000., -10., 10.);
    TH2F* h_delta_z_vs_t_corr = new TH2F("h_delta_z_vs_t_corr", "|#DeltaZ| vs t_{#gamma}^{corrected}", 6, 0, 6, 6, 0, 6);
    TH2F* h_delta_z_vs_t_corr_loose = new TH2F("h_delta_z_vs_t_corr_loose", "|#DeltaZ| vs t_{#gamma}^{corrected}", 6, 0, 6, 6, 0, 6);

    TProfile* p_t_vs_et_fit = new TProfile("p_t_vs_e_fit", "Profile of t_{#gamma}^{corrected} vs E, Medium Gain, EMB, Loose Photons", 100, 0., 2000., -10., 10.);
    TProfile* p_t_vs_et_fit_val = new TProfile("p_t_vs_e_fit_val", "Profile of t_{#gamma}^{corrected} vs E, Medium Gain, EMB, Loose Photons", 100, 0., 2000., -10., 10.);
    TProfile* p_t_vs_et_tightPh = new TProfile("p_t_vs_e_tightPh", "Profile of t_{#gamma}^{corrected} vs E, Medium Gain, EMB, Loose Photons", 100, 0., 2000., -10., 10.);

    TProfile* p_t_vs_et_fit_tightPh = new TProfile("p_t_vs_e_fit_tightPh", "Profile of t_{#gamma}^{corrected} vs E, Medium Gain, EMB, Loose Photons", 100, 0., 2000., -10., 10.);
    TProfile* p_t_vs_et_fit_val_tightPh = new TProfile("p_t_vs_e_fit_val_tightPh", "Profile of t_{#gamma}^{corrected} vs E, Medium Gain, EMB, Loose Photons", 100, 0., 2000., -10., 10.);

    TH1F* h_ph_t_uncor_val = new TH1F("h_ph_t_uncor_val", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_t_uncor_tightPh = new TH1F("h_ph_t_uncor_tightPh", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_t_uncor_val_tightPh = new TH1F("h_ph_t_uncor_val_tightPh", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);

    //General
    TH1F* h_ph_e = new TH1F("h_ph_e", "Leading photon energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_et = new TH1F("h_ph_et", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta = new TH1F("h_ph_eta", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi = new TH1F("h_ph_phi", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor = new TH1F("h_ph_t_uncor", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_t_corr = new TH1F("h_ph_t_corr", "t_{#gamma}^{(corrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_gain = new TH1F("h_ph_gain", "Gain", 3, 0, 3);
    TH1F* h_ph_delta_z = new TH1F("h_ph_delta_z", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_ph_dphi = new TH1F("h_ph_dphi", "#Delta#phi", 100, -1., 1.);
    TH1F* h_MET = new TH1F("h_MET", "E_{T}^{Miss} [GeV]", 1000, 0., 1000.);
    TH1F* h_MET_tight = new TH1F("h_MET_tight", "E_{T}^{Miss} [GeV]", 1000, 0., 1000.);
    TH2F* h_t_vs_et = new TH2F("h_t_vs_et", "t_{#gamma} vs E_{T}", 100, 0., 2000., 100, -10., 10.);

    //EMB
    TH1F* h_ph_et_emb = new TH1F("h_ph_et_emb", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_emb = new TH1F("h_ph_eta_emb", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emb = new TH1F("h_ph_phi_emb", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emb = new TH1F("h_ph_t_uncor_emb", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_gain_emb = new TH1F("h_ph_gain_emb", "Gain", 3, 0, 3);
    TH1F* h_ph_delta_z_emb = new TH1F("h_ph_delta_z_emb", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emb = new TH1F("h_MET_emb", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //EMEC
    TH1F* h_ph_et_emec = new TH1F("h_ph_et_emec", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_emec = new TH1F("h_ph_eta_emec", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emec = new TH1F("h_ph_phi_emec", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emec = new TH1F("h_ph_t_uncor_emec", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_gain_emec = new TH1F("h_ph_gain_emec", "Gain", 3, 0, 3);
    TH1F* h_ph_delta_z_emec = new TH1F("h_ph_delta_z_emec", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emec = new TH1F("h_MET_emec", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //HI-GAIN
    TH1F* h_ph_et_higain = new TH1F("h_ph_et_higain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_higain = new TH1F("h_ph_eta_higain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_higain = new TH1F("h_ph_phi_higain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_higain = new TH1F("h_ph_t_uncor_higain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_maxEcell_E_higain = new TH1F("h_ph_maxEcell_E_higain", "Max Ecell Energy [GeV]", 1000, 0., 300.);
    TH1F* h_ph_delta_z_higain = new TH1F("h_ph_delta_z_higain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_higain = new TH1F("h_MET_higain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //MED-GAIN
    TH1F* h_ph_et_medgain = new TH1F("h_ph_et_medgain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_medgain = new TH1F("h_ph_eta_medgain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_medgain = new TH1F("h_ph_phi_medgain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_medgain = new TH1F("h_ph_t_uncor_medgain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_maxEcell_E_medgain = new TH1F("h_ph_maxEcell_E_medgain", "Max Ecell Energy [GeV]", 1000, 0., 500.);
    TH1F* h_ph_delta_z_medgain = new TH1F("h_ph_delta_z_medgain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_medgain = new TH1F("h_MET_medgain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //LOW-GAIN
    TH1F* h_ph_et_logain = new TH1F("h_ph_et_logain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_logain = new TH1F("h_ph_eta_logain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_logain = new TH1F("h_ph_phi_logain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_logain = new TH1F("h_ph_t_uncor_logain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_maxEcell_E_logain = new TH1F("h_ph_maxEcell_E_logain", "Max Ecell Energy [GeV]", 1000, 0., 1000.);
    TH1F* h_ph_delta_z_logain = new TH1F("h_ph_delta_z_logain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_logain = new TH1F("h_MET_logain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //HI-GAIN+EMB
    TH1F* h_ph_et_emb_higain = new TH1F("h_ph_et_emb_higain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_emb_higain = new TH1F("h_ph_eta_emb_higain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emb_higain = new TH1F("h_ph_phi_emb_higain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emb_higain = new TH1F("h_ph_t_uncor_emb_higain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_delta_z_emb_higain = new TH1F("h_ph_delta_z_emb_higain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emb_higain = new TH1F("h_MET_emb_higain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //HI-GAIN+EMEC
    TH1F* h_ph_et_emec_higain = new TH1F("h_ph_et_emec_higain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_emec_higain = new TH1F("h_ph_eta_emec_higain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emec_higain = new TH1F("h_ph_phi_emec_higain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emec_higain = new TH1F("h_ph_t_uncor_emec_higain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_delta_z_emec_higain = new TH1F("h_ph_delta_z_emec_higain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emec_higain = new TH1F("h_MET_emec_higain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //MED-GAIN+EMB
    TH1F* h_ph_et_emb_medgain = new TH1F("h_ph_et_emb_medgain", "Leading photon transverse energy [GeV]", 100, 0., 1500.);
    TH1F* h_ph_eta_emb_medgain = new TH1F("h_ph_eta_emb_medgain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emb_medgain = new TH1F("h_ph_phi_emb_medgain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emb_medgain = new TH1F("h_ph_t_uncor_emb_medgain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_delta_z_emb_medgain = new TH1F("h_ph_delta_z_emb_medgain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emb_medgain = new TH1F("h_MET_emb_medgain", "E_{T}^{Miss} [GeV]", 100, 0., 1000.);
    TH1F* h_MET_emb_medgain_iso = new TH1F("h_MET_emb_medgain_iso", "E_{T}^{Miss} [GeV]", 500, 0., 1000.);

    //MED-GAIN+EMEC
    TH1F* h_ph_et_emec_medgain = new TH1F("h_ph_et_emec_medgain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_emec_medgain = new TH1F("h_ph_eta_emec_medgain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emec_medgain = new TH1F("h_ph_phi_emec_medgain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emec_medgain = new TH1F("h_ph_t_uncor_emec_medgain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_delta_z_emec_medgain = new TH1F("h_ph_delta_z_emec_medgain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emec_medgain = new TH1F("h_MET_emec_medgain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //LOW-GAIN+EMB
    TH1F* h_ph_et_emb_logain = new TH1F("h_ph_et_emb_logain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_emb_logain = new TH1F("h_ph_eta_emb_logain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emb_logain = new TH1F("h_ph_phi_emb_logain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emb_logain = new TH1F("h_ph_t_uncor_emb_logain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_delta_z_emb_logain = new TH1F("h_ph_delta_z_emb_logain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emb_logain = new TH1F("h_MET_emb_logain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    //LOW-GAIN+EMEC
    TH1F* h_ph_et_emec_logain = new TH1F("h_ph_et_emec_logain", "Leading photon transverse energy [GeV]", 1000, 0., 1500.);
    TH1F* h_ph_eta_emec_logain = new TH1F("h_ph_eta_emec_logain", "Leading photon #eta", 1000, -3., 3.);
    TH1F* h_ph_phi_emec_logain = new TH1F("h_ph_phi_emec_logain", "Leading photon #phi", 1000, -3.15, 3.15);
    TH1F* h_ph_t_uncor_emec_logain = new TH1F("h_ph_t_uncor_emec_logain", "t_{#gamma}^{(uncorrected)} [ns]", 1000, -10., 10.);
    TH1F* h_ph_delta_z_emec_logain = new TH1F("h_ph_delta_z_emec_logain", "#Delta Z [mm]", 1000, -1000., 1000.);
    TH1F* h_MET_emec_logain = new TH1F("h_MET_emec_logain", "E_{T}^{Miss} [GeV]", 1000, 0., 250.);

    double delta_z;
    double dphi;
    double met_phi;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        if (jentry % 100000 == 0) cout << "Processed " << jentry << " events" << endl;

        //START FILL PLOTS

        if(fabs(ph_eta1) < 1.37) emb = true; //BARREL
        if((fabs(ph_eta1) > 1.52) && (fabs(ph_eta1) < 2.5)) emec = true; //ENDCAP
        if(ph1_maxEcell_gain == 0) higain = true; //HIGAIN
        if(ph1_maxEcell_gain == 1) medgain = true; //MEDGAIN
        if(ph1_maxEcell_gain == 2) logain = true; //LOWGAIN
        if( (ph1_topoetcone20/1000. < 0.065*ph_pt1) && (ph1_ptcone20/(ph_pt1*1000.) < 0.05) ) looseIso = true; //Loose Iso

        nTot++;

        h_MET->Fill(m_met/1000.);
        if((n_ph == 1) && (ph_pt1 > 150.) && emb ){

            met_phi = TMath::ATan2(m_mpy, m_mpx);

            if(fabs(ph_phi1 - met_phi) > TMath::Pi()) dphi = 2*TMath::Pi() - fabs(ph_phi1 - met_phi);
            else dphi = fabs(ph_phi1 - met_phi);

            h_ph_dphi->Fill(TMath::Cos(dphi));

            delta_z = ph_calo_z1 - PV_z;
            if( (-4. <= ph_t1) && (ph_t1 < 0.5) ){
                if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr->Fill(0., (int)(fabs(delta_z)/40));
                else h_delta_z_vs_t_corr->Fill(0., 5); 
            }
            if( (0.5 <= ph_t1) && (ph_t1 < 1.1) ){
                if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr->Fill(1, (int)(fabs(delta_z)/40));
                else h_delta_z_vs_t_corr->Fill(1, 5); 
            }
            if( (1.1 <= ph_t1) && (ph_t1 < 1.3) ){
                if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr->Fill(2, (int)(fabs(delta_z)/40));
                else h_delta_z_vs_t_corr->Fill(2, 5); 
            }
            if( (1.3 <= ph_t1) && (ph_t1 < 1.5) ){
                if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr->Fill(3, (int)(fabs(delta_z)/40));
                else h_delta_z_vs_t_corr->Fill(3, 5); 
            }
            if( (1.5 <= ph_t1) && (ph_t1 < 1.8) ){
                if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr->Fill(4, (int)(fabs(delta_z)/40));
                else h_delta_z_vs_t_corr->Fill(4, 5); 
            }
            if( (1.8 <= ph_t1) && (ph_t1 < 4.) ){
                if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr->Fill(5, (int)(fabs(delta_z)/40));
                else h_delta_z_vs_t_corr->Fill(5, 5); 
            }

            if(looseIso){
                delta_z = ph_calo_z1 - PV_z;
                if( (-4. <= ph_t1) && (ph_t1 < 0.5) ){
                    if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr_loose->Fill(0., (int)(fabs(delta_z)/40));
                    else h_delta_z_vs_t_corr_loose->Fill(0., 5); 
                }
                if( (0.5 <= ph_t1) && (ph_t1 < 1.1) ){
                    if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr_loose->Fill(1, (int)(fabs(delta_z)/40));
                    else h_delta_z_vs_t_corr_loose->Fill(1, 5); 
                }
                if( (1.1 <= ph_t1) && (ph_t1 < 1.3) ){
                    if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr_loose->Fill(2, (int)(fabs(delta_z)/40));
                    else h_delta_z_vs_t_corr_loose->Fill(2, 5); 
                }
                if( (1.3 <= ph_t1) && (ph_t1 < 1.5) ){
                    if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr_loose->Fill(3, (int)(fabs(delta_z)/40));
                    else h_delta_z_vs_t_corr_loose->Fill(3, 5); 
                }
                if( (1.5 <= ph_t1) && (ph_t1 < 1.8) ){
                    if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr_loose->Fill(4, (int)(fabs(delta_z)/40));
                    else h_delta_z_vs_t_corr_loose->Fill(4, 5); 
                }
                if( (1.8 <= ph_t1) && (ph_t1 < 4.) ){
                    if(fabs(delta_z) < 200.) h_delta_z_vs_t_corr_loose->Fill(5, (int)(fabs(delta_z)/40));
                    else h_delta_z_vs_t_corr_loose->Fill(5, 5); 
                }
            } 

            h_ph_t_uncor->Fill(ph1_maxEcell_time);
            h_ph_t_corr->Fill(ph_t1);
            h_ph_et->Fill(ph_pt1); 
            h_ph_eta->Fill(ph_eta1); 
            h_ph_phi->Fill(ph_phi1); 

            h_MET->Fill(m_met/1000.);
            h_ph_delta_z->Fill(ph_calo_z1 - PV_z);

            if(higain){
                h_ph_t_uncor_higain->Fill(ph1_maxEcell_time);
                h_ph_et_higain->Fill(ph_pt1); 
                h_ph_eta_higain->Fill(ph_eta1); 
                h_ph_phi_higain->Fill(ph_phi1); 

                h_ph_maxEcell_E_higain->Fill(ph1_maxEcell_E);
                h_MET_higain->Fill(m_met/1000.);
                h_ph_delta_z_higain->Fill(ph_calo_z1 - PV_z);
                if(emb){
                    h_ph_et_emb_higain->Fill(ph_pt1); 
                    h_ph_eta_emb_higain->Fill(ph_eta1); 
                    h_ph_phi_emb_higain->Fill(ph_phi1); 

                    h_ph_t_uncor_emb_higain->Fill(ph1_maxEcell_time);
                    h_MET_emb_higain->Fill(m_met/1000.);
                    h_ph_delta_z_emb_higain->Fill(ph_calo_z1 - PV_z);
                }
                if(emec){
                    h_ph_et_emec_higain->Fill(ph_pt1); 
                    h_ph_eta_emec_higain->Fill(ph_eta1); 
                    h_ph_phi_emec_higain->Fill(ph_phi1); 

                    h_ph_t_uncor_emec_higain->Fill(ph1_maxEcell_time);
                    h_MET_emec_higain->Fill(m_met/1000.);
                    h_ph_delta_z_emec_higain->Fill(ph_calo_z1 - PV_z);
                }
            }
            if(medgain){
                h_ph_t_uncor_medgain->Fill(ph1_maxEcell_time);
                h_ph_et_medgain->Fill(ph_pt1); 
                h_ph_eta_medgain->Fill(ph_eta1); 
                h_ph_phi_medgain->Fill(ph_phi1); 

                h_ph_maxEcell_E_medgain->Fill(ph1_maxEcell_E);
                h_MET_medgain->Fill(m_met/1000.);
                h_ph_delta_z_medgain->Fill(ph_calo_z1 - PV_z);
                if(emb){
                    h_ph_et_emb_medgain->Fill(ph_pt1); 
                    h_ph_eta_emb_medgain->Fill(ph_eta1); 
                    h_ph_phi_emb_medgain->Fill(ph_phi1); 

                    h_ph_t_uncor_emb_medgain->Fill(ph1_maxEcell_time);
                    h_MET_emb_medgain->Fill(m_met/1000.);
                    h_ph_delta_z_emb_medgain->Fill(ph_calo_z1 - PV_z);
                }
                if(emec){
                    h_ph_et_emec_medgain->Fill(ph_pt1); 
                    h_ph_eta_emec_medgain->Fill(ph_eta1); 
                    h_ph_phi_emec_medgain->Fill(ph_phi1); 

                    h_ph_t_uncor_emec_medgain->Fill(ph1_maxEcell_time);
                    h_MET_emec_medgain->Fill(m_met/1000.);
                    h_ph_delta_z_emec_medgain->Fill(ph_calo_z1 - PV_z);
                }
            }
            if(logain){
                h_ph_t_uncor_logain->Fill(ph1_maxEcell_time);
                h_ph_et_logain->Fill(ph_pt1); 
                h_ph_eta_logain->Fill(ph_eta1); 
                h_ph_phi_logain->Fill(ph_phi1); 

                h_ph_maxEcell_E_logain->Fill(ph1_maxEcell_E);
                h_MET_logain->Fill(m_met/1000.);
                h_ph_delta_z_logain->Fill(ph_calo_z1 - PV_z);

                if(emb){
                    h_ph_et_emb_logain->Fill(ph_pt1); 
                    h_ph_eta_emb_logain->Fill(ph_eta1); 
                    h_ph_phi_emb_logain->Fill(ph_phi1); 

                    h_ph_t_uncor_emb_logain->Fill(ph1_maxEcell_time);
                    h_MET_emb_logain->Fill(m_met/1000.);
                    h_ph_delta_z_emb_logain->Fill(ph_calo_z1 - PV_z);
                }
                if(emec){
                    h_ph_et_emec_logain->Fill(ph_pt1); 
                    h_ph_eta_emec_logain->Fill(ph_eta1); 
                    h_ph_phi_emec_logain->Fill(ph_phi1); 

                    h_ph_t_uncor_emec_logain->Fill(ph1_maxEcell_time);
                    h_MET_emec_logain->Fill(m_met/1000.);
                    h_ph_delta_z_emec_logain->Fill(ph_calo_z1 - PV_z);
                }
            }
        }//END FILL PLOTS


        emb = false;
        emec = false;
        logain = false;
        medgain = false;
        higain = false;
    }//END LOOP

    TFile f(Form("NtuplePlots/%s_plots.root", fileName.c_str()),"recreate");

    h_delta_z_vs_t_corr->Write();
    h_delta_z_vs_t_corr_loose->Write();

    h_ph_e->Write();
    h_ph_et->Write();
    h_ph_eta->Write();
    h_ph_phi->Write();
    h_ph_t_uncor->Write();
    h_ph_t_corr->Write();
    h_ph_gain->Write();
    h_ph_delta_z->Write();
    h_ph_dphi->Write();
    h_MET->Write();
    h_MET_tight->Write();
    h_t_vs_et->Write();

    //EMB
    h_ph_et_emb->Write();
    h_ph_eta_emb->Write();
    h_ph_phi_emb->Write();
    h_ph_t_uncor_emb->Write();
    h_ph_gain_emb->Write();
    h_ph_delta_z_emb->Write();
    h_MET_emb->Write();

    //EMEC
    h_ph_et_emec->Write();
    h_ph_eta_emec->Write();
    h_ph_phi_emec->Write();
    h_ph_t_uncor_emec->Write();
    h_ph_gain_emec->Write();
    h_ph_delta_z_emec->Write();
    h_MET_emec->Write();

    //HIGAIN
    h_ph_et_higain->Write();
    h_ph_eta_higain->Write();
    h_ph_phi_higain->Write();
    h_ph_t_uncor_higain->Write();
    h_ph_maxEcell_E_higain->Write();
    h_ph_delta_z_higain->Write();
    h_MET_higain->Write();

    //MEDGAIN
    h_ph_et_medgain->Write();
    h_ph_eta_medgain->Write();
    h_ph_phi_medgain->Write();
    h_ph_t_uncor_medgain->Write();
    h_ph_maxEcell_E_medgain->Write();
    h_ph_delta_z_medgain->Write();
    h_MET_medgain->Write();

    //LOGAIN
    h_ph_et_logain->Write();
    h_ph_eta_logain->Write();
    h_ph_phi_logain->Write();
    h_ph_t_uncor_logain->Write();
    h_ph_maxEcell_E_logain->Write();
    h_ph_delta_z_logain->Write();
    h_MET_logain->Write();

    //HIGAIN+EMB
    h_ph_et_emb_higain->Write();
    h_ph_eta_emb_higain->Write();
    h_ph_phi_emb_higain->Write();
    h_ph_t_uncor_emb_higain->Write();
    h_ph_delta_z_emb_higain->Write();
    h_MET_emb_higain->Write();

    //HIGAIN+EMEC
    h_ph_et_emec_higain->Write();
    h_ph_eta_emec_higain->Write();
    h_ph_phi_emec_higain->Write();
    h_ph_t_uncor_emec_higain->Write();
    h_ph_delta_z_emec_higain->Write();
    h_MET_emec_higain->Write();

    //MEDGAIN+EMB
    h_ph_et_emb_medgain->Write();
    h_ph_eta_emb_medgain->Write();
    h_ph_phi_emb_medgain->Write();
    h_ph_t_uncor_emb_medgain->Write();
    h_ph_delta_z_emb_medgain->Write();
    h_MET_emb_medgain->Write();
    h_MET_emb_medgain_iso->Write();

    //MEDGAIN+EMEC
    h_ph_et_emec_medgain->Write();
    h_ph_eta_emec_medgain->Write();
    h_ph_phi_emec_medgain->Write();
    h_ph_t_uncor_emec_medgain->Write();
    h_ph_delta_z_emec_medgain->Write();
    h_MET_emec_medgain->Write();

    //LOGAIN+EMB
    h_ph_et_emb_logain->Write();
    h_ph_eta_emb_logain->Write();
    h_ph_phi_emb_logain->Write();
    h_ph_t_uncor_emb_logain->Write();
    h_ph_delta_z_emb_logain->Write();
    h_MET_emb_logain->Write();

    //LOGAIN+EMEC
    h_ph_et_emec_logain->Write();
    h_ph_eta_emec_logain->Write();
    h_ph_phi_emec_logain->Write();
    h_ph_t_uncor_emec_logain->Write();
    h_ph_delta_z_emec_logain->Write();
    h_MET_emec_logain->Write();

    f.Close();
}
