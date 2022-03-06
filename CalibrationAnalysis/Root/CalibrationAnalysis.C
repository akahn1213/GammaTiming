#define CalibrationAnalysis_cxx
#include "../CalibrationAnalysis/CalibrationAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



void CalibrationAnalysis(){
    class CalibrationAnalysis t;
    t.Loop();
}

void CalibrationAnalysis::Loop()
{
    //   In a ROOT session, you can do:
    //      root> .L CalibrationAnalysis.C
    //      root> CalibrationAnalysis t
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

    //HISTOGRAMS
    bookHists();


    double t_pass0;
    double t_pass1;
    double t_pass2;
    double t_pass3;
    double t_pass4;
    int febN;
    int slot = -1;
    unsigned long int onlId;
    double e;

    string line;

    TFile* corrections_file = TFile::Open("../Data/corrections/2016_corrections.root");


    //Read inputs
    //int maxentries = 500000;
    int maxentries = -1;

    Long64_t nentries = fChain->GetEntriesFast();
    if(maxentries < 0) maxentries = nentries;
    Long64_t nbytes = 0, nb = 0;
    //PASS 0 LOOP
    for (Long64_t jentry=0; jentry<maxentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( (jentry > 0) && (jentry%100000 == 0) ) cout << "Pass 0: Events Processed: " << jentry << endl;
        if(badOnlidCheck(ph1_maxEcell_onlId)) continue;

        if((ph1_maxEcell_gain == 2) && (m_met < 100.) && (ph_pt1 > 150.) && (n_ph == 1)){
            h_t_uncorr->Fill(ph1_maxEcell_time);
            t_pass0 = tofCorrection(ph1_maxEcell_x, ph1_maxEcell_y, ph1_maxEcell_z, ph1_maxEcell_time, PV_z); 

            if(fabs(ph_eta1) < 1.37){
                if( ph_eta1 < 0){//EMBA
                    h_t_pass0_emba->Fill(t_pass0);
                    p_t_pass0_vs_runN_emba->Fill(runNumber, t_pass0);
                }
                if(ph_eta1 > 0){//EMBC
                    h_t_pass0_embc->Fill(t_pass0);
                    p_t_pass0_vs_runN_embc->Fill(runNumber, t_pass0);
                }
                h_t_pass0_emb->Fill(t_pass0);
            }//BARREL

            if((fabs(ph_eta1) > 1.52) && (fabs(ph_eta1) < 2.5)){
                if( ph_eta1 < 0){//EMECA
                    h_t_pass0_emeca->Fill(t_pass0);
                    p_t_pass0_vs_runN_emeca->Fill(runNumber, t_pass0);
                }
                if(ph_eta1 > 0){//EMECC
                    h_t_pass0_emecc->Fill(t_pass0);
                    p_t_pass0_vs_runN_emecc->Fill(runNumber, t_pass0);
                }

                h_t_pass0_emec->Fill(t_pass0);
            }//ENDCAP

            h_t_pass0->Fill(t_pass0);
            p_t_pass0_vs_runN->Fill(runNumber, t_pass0);

        }
    }//END PASS 0 LOOP
    nbytes = 0;
    nb = 0;

    //PASS 1 LOOP
    for (Long64_t jentry=0; jentry<maxentries;jentry++) {//EVENT LOOP
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( (jentry > 0) && (jentry%100000 == 0) ) cout << "Pass 1: Events Processed: " << jentry << endl;
        if(badOnlidCheck(ph1_maxEcell_onlId)) continue;

        if((ph1_maxEcell_gain == 2) && (m_met < 100.) && (ph_pt1 > 150.) && (n_ph == 1)){

            t_pass0 = tofCorrection(ph1_maxEcell_x, ph1_maxEcell_y, ph1_maxEcell_z, ph1_maxEcell_time, PV_z); 

            if(fabs(ph_eta1) < 1.37){
                if( ph_eta1 < 0){//EMBA
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emba->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                    h_t_pass1_emba->Fill(t_pass1);
                    p_t_pass1_vs_e_emba->Fill(ph_pt1*TMath::CosH(ph_eta1), t_pass1);
                    p_t_pass1_vs_runN_emba->Fill(runNumber, t_pass1);
                }
                if(ph_eta1 > 0){//EMBC
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_embc->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                    h_t_pass1_embc->Fill(t_pass1);
                    p_t_pass1_vs_e_embc->Fill(ph_pt1*TMath::CosH(ph_eta1), t_pass1);
                    p_t_pass1_vs_runN_embc->Fill(runNumber, t_pass1);
                }

                h_t_pass1_emb->Fill(t_pass1);
            }//BARREL

            if((fabs(ph_eta1) > 1.52) && (fabs(ph_eta1) < 2.5)){
                if( ph_eta1 < 0){//EMECA
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emeca->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                    h_t_pass1_emeca->Fill(t_pass1);
                    p_t_pass1_vs_e_emeca->Fill(ph_pt1*TMath::CosH(ph_eta1), t_pass1);
                    p_t_pass1_vs_runN_emeca->Fill(runNumber, t_pass1);
                }
                if(ph_eta1 > 0){//EMECC
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emecc->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                    h_t_pass1_emecc->Fill(t_pass1);
                    p_t_pass1_vs_e_emecc->Fill(ph_pt1*TMath::CosH(ph_eta1), t_pass1);
                    p_t_pass1_vs_runN_emecc->Fill(runNumber, t_pass1);
                }

                h_t_pass1_emec->Fill(t_pass1);
            }//ENDCAP

            onlId = ph1_maxEcell_onlId>>32; 
            febN = getFebNo(onlId);

            p_t_pass1_vs_febN->Fill(febN, t_pass1);

            h_t_pass1->Fill(t_pass1);
            p_t_pass1_vs_runN->Fill(runNumber, t_pass1);

        }
    }//END PASS 1 LOOP

    nbytes = 0;
    nb = 0;
    //PASS 2 LOOP
    for (Long64_t jentry=0; jentry<maxentries;jentry++) {//EVENT LOOP
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( (jentry > 0) && (jentry%100000 == 0) ) cout << "Pass 2: Events Processed: " << jentry << endl;
        if(badOnlidCheck(ph1_maxEcell_onlId)) continue;

        if((ph1_maxEcell_gain == 2) && (m_met < 100.) && (ph_pt1 > 150.) && (n_ph == 1)){

            t_pass0 = tofCorrection(ph1_maxEcell_x, ph1_maxEcell_y, ph1_maxEcell_z, ph1_maxEcell_time, PV_z); 

            if(fabs(ph_eta1) < 1.37){
                if( ph_eta1 < 0){//EMBA
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emba->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMBC
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_embc->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }

            }//BARREL

            if((fabs(ph_eta1) > 1.52) && (fabs(ph_eta1) < 2.5)){
                if( ph_eta1 < 0){//EMECA
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emeca->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMECC
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emecc->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
            }//ENDCAP

            onlId = ph1_maxEcell_onlId>>32; 
            febN = getFebNo(onlId);
            t_pass2 = t_pass1 - (p_t_pass1_vs_febN->GetBinContent(febN));
            p_t_pass2_vs_febN->Fill(febN, t_pass2);
            h_t_pass2->Fill(t_pass2);
            if( (febN < 62) || ( (febN > 130) && (febN < 190) ) ) h_t_pass2_narrow->Fill(t_pass2);
            if(  ( (febN > 64) && (febN < 128) ) || ( (febN > 195) && (febN < 256) ) ) h_t_pass2_wide->Fill(t_pass2);

            h_e_fullFEB->Fill(ph_pt1*TMath::CosH(ph_eta1));

            for(int i = 1; i < 257; i++){
                if(febN == i){
                    h_t_fn[i-1]->Fill(t_pass1);
                    if(  (((int)(i/64)) % 2) == 0 ){//Narrow FEB Plots            
                        h_t_febs_narrow->Fill(t_pass1);
                        h_t_corr_febs_narrow->Fill(t_pass2);
                        h_e_narrow->Fill(ph_pt1*TMath::CosH(ph_eta1));
                        h_e_cell_narrow->Fill(ph1_maxEcell_E/1000);
                        h_et_narrow->Fill(ph_pt1);
                        h_eta_narrow->Fill(ph_eta1);
                        h_phi_narrow->Fill(ph_phi1);
                    }
                    else{//Wide FEB Plots
                        h_t_febs_wide->Fill(t_pass1);
                        h_t_corr_febs_wide->Fill(t_pass2);
                        h_e_wide->Fill(ph_pt1*TMath::CosH(ph_eta1));
                        h_e_cell_wide->Fill(ph1_maxEcell_E/1000);
                        h_et_wide->Fill(ph_pt1);
                        h_eta_wide->Fill(ph_eta1);
                        h_phi_wide->Fill(ph_phi1);
                    }
                }
            }


            slot = (int)(febN/32);
            if(slot < 8){
                h_t_febs_slot[slot]->Fill(t_pass1);
                h_t_corr_febs_slot[slot]->Fill(t_pass2);
                h_e_slot[slot]->Fill(ph_pt1*TMath::CosH(ph_eta1));
                h_e_cell_slot[slot]->Fill(ph1_maxEcell_E/1000);
                h_et_slot[slot]->Fill(ph_pt1);
                h_eta_slot[slot]->Fill(ph_eta1);
                h_phi_slot[slot]->Fill(ph_phi1);
                p_t_vs_E_slot[slot]->Fill(ph_pt1*TMath::CosH(ph_eta1), t_pass1);
                p_t_corr_vs_E_slot[slot]->Fill(ph_pt1*TMath::CosH(ph_eta1), t_pass2);
                p_t_vs_ET_slot[slot]->Fill(ph_pt1, t_pass1);
                p_t_corr_vs_ET_slot[slot]->Fill(ph_pt1, t_pass2);
                p_t_vs_E_cell_slot[slot]->Fill(ph1_maxEcell_E/1000, t_pass1);
                p_t_pass2_vs_E_cell_slot[slot]->Fill(ph1_maxEcell_E/1000, t_pass2);

            }



            if(febN > 256){
                h_t_febs_emec->Fill(t_pass1);
                h_t_corr_febs_emec->Fill(t_pass2);
                h_e_emec->Fill(ph_pt1*TMath::CosH(ph_eta1));
                h_e_cell_emec->Fill(ph1_maxEcell_E/1000);
                h_et_emec->Fill(ph_pt1);
                h_eta_emec->Fill(ph_eta1);
                h_phi_emec->Fill(ph_phi1);
                p_t_corr_emec_vs_subslot->Fill(getSubslotNo(onlId), t_pass2);
            }
        }
    }//END PASS 2 LOOP


    nbytes = 0;
    nb = 0;
    //PASS 3 LOOP
    for (Long64_t jentry=0; jentry<maxentries;jentry++) {//EVENT LOOP
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( (jentry > 0) && (jentry%100000 == 0) ) cout << "Pass 3: Events Processed: " << jentry << endl;
        if(badOnlidCheck(ph1_maxEcell_onlId)) continue;

        if((ph1_maxEcell_gain == 2) && (m_met < 100.) && (ph_pt1 > 150.) && (n_ph == 1)){

            t_pass0 = tofCorrection(ph1_maxEcell_x, ph1_maxEcell_y, ph1_maxEcell_z, ph1_maxEcell_time, PV_z); 

            if(fabs(ph_eta1) < 1.37){
                if( ph_eta1 < 0){//EMBA
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emba->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMBC
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_embc->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }

            }//BARREL

            if((fabs(ph_eta1) > 1.52) && (fabs(ph_eta1) < 2.5)){
                if( ph_eta1 < 0){//EMECA
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emeca->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMECC
                    t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emecc->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }

            }//ENDCAP

            onlId = ph1_maxEcell_onlId>>32; 
            febN = getFebNo(onlId);
            t_pass2 = t_pass1 - (p_t_pass1_vs_febN->GetBinContent(febN));

            e = ph1_maxEcell_E/1000;
            slot = (int)(febN/32);
            t_pass3 = -999.;
            switch(slot){
                case 0: {
                            if(e < 650.) t_pass3 = t_pass2 - (1.34131e-06*e*e - 2.59592e-03*e + 5.98932e-01);                       
                            else t_pass3 = t_pass2 - (1.34131e-06*650*650 - 2.59592e-03*650 + 5.98932e-01);                       
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[0]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
                case 1: {
                            if(e < 650.) t_pass3 = t_pass2 - (2.30630e-06*e*e - 3.59088e-03*e + 8.37710e-01);                       
                            else t_pass3 = t_pass2 - (2.30630e-06*650*650 - 3.59088e-03*650 + 8.37710e-01);                       
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[1]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
                case 2: {
                            if(e < 700.) t_pass3 = t_pass2 - (1.09187e-06*e*e - 1.85810e-03*e + 6.16434e-01);                       
                            else t_pass3 = t_pass2 - (1.09187e-06*700*700 - 1.85810e-03*700 + 6.16434e-01);                       
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[2]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
                case 3: {
                            if(e < 700.) t_pass3 = t_pass2 - (-5.14276e-06*e*e + 4.32677e-03*e - 8.73504e-01);                        
                            else t_pass3 = t_pass2 - (-5.14276e-06*700*700 + 4.32677e-03*700 - 8.73504e-01);                        
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[3]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
                case 4: {
                            if(e < 650.) t_pass3 = t_pass2 - (1.93255e-06*e*e - 3.17103e-03*e + 7.11106e-01);                       
                            else t_pass3 = t_pass2 - (1.93255e-06*650*650 - 3.17103e-03*650 + 7.11106e-01);                       
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[0]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
                case 5: {
                            if(e < 650.) t_pass3 = t_pass2 - (-2.52727e-07*e*e - 1.66121e-03*e + 4.97491e-01);                        
                            else t_pass3 = t_pass2 - (-2.52727e-07*650*650 - 1.66121e-03*650 + 4.97491e-01);                        
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[1]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
                case 6: {
                            if(e < 700.) t_pass3 = t_pass2 - (4.57446e-06*e*e - 5.24988e-03*e + 1.41012e+00);                       
                            else t_pass3 = t_pass2 - (4.57446e-06*700*700 - 5.24988e-03*700 + 1.41012e+00);                       
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[2]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
                case 7: {
                            if(e < 700.) t_pass3 = t_pass2 - (-5.90292e-06*e*e + 5.10021e-03*e - 1.06100e+00);                        
                            else t_pass3 = t_pass2 - (-5.90292e-06*700*700 + 5.10021e-03*700 - 1.06100e+00);                        
                            h_t_pass3_slot[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot[3]->Fill(getChannelNo(onlId)%128, t_pass3);
                            break;
                        }
            }
            if(t_pass3 > -999.) h_t_pass3->Fill(t_pass3);
            if(fabs(ph_eta1) < 1.37) h_t_pass3_emb->Fill(t_pass3);
        }
    }//END PASS 3 LOOP

    nbytes = 0;
    nb = 0;
    //PASS 4 LOOP
    for (Long64_t jentry=0; jentry<maxentries;jentry++) {//EVENT LOOP
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        if( (jentry > 0) && (jentry%100000 == 0) ) cout << "Pass 4: Events Processed: " << jentry << endl;
        if(badOnlidCheck(ph1_maxEcell_onlId)) continue;

        if((ph1_maxEcell_gain == 2) && (m_met < 100.) && (ph_pt1 > 150.) && (n_ph == 1)){

            //APPLY PASS 0 CORRECTIONS
            t_pass0 = tofCorrection(ph1_maxEcell_x, ph1_maxEcell_y, ph1_maxEcell_z, ph1_maxEcell_time, PV_z); 

            //APPLY PASS 1 CORRECTIONS
            if(fabs(ph_eta1) < 1.37){
                if( ph_eta1 < 0){//EMBA
                    t_pass1 = t_pass0 - (((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_emba"))->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMBC
                    t_pass1 = t_pass0 - (((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_embc"))->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }

            }//BARREL

            if((fabs(ph_eta1) > 1.52) && (fabs(ph_eta1) < 2.5)){
                if( ph_eta1 < 0){//EMECA
                    t_pass1 = t_pass0 - (((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_emeca"))->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMECC
                    t_pass1 = t_pass0 - (((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_emecc"))->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }

            }//ENDCAP

            //APPLY PASS 2 CORRECTIONS
            onlId = ph1_maxEcell_onlId>>32; 
            febN = getFebNo(onlId);
            t_pass2 = t_pass1 - (((TProfile*)corrections_file->Get("p_t_pass1_vs_febN"))->GetBinContent(febN));

            //APPLY PASS 3 & 4 CORRECTIONS  
            e = ph1_maxEcell_E/1000;
            slot = (int)(febN/32);
            t_pass3 = -999.;
            t_pass4 = -999.;
            switch(slot){
                case 0: {
                            if(e < 650.) t_pass3 = t_pass2 - (1.34131e-06*e*e - 2.59592e-03*e + 5.98932e-01);                       
                            else t_pass3 = t_pass2 - (1.34131e-06*650*650 - 2.59592e-03*650 + 5.98932e-01);                       

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
                case 1: {
                            if(e < 650.) t_pass3 = t_pass2 - (2.30630e-06*e*e - 3.59088e-03*e + 8.37710e-01);                       
                            else t_pass3 = t_pass2 - (2.30630e-06*650*650 - 3.59088e-03*650 + 8.37710e-01);                       

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
                case 2: {
                            if(e < 700.) t_pass3 = t_pass2 - (1.09187e-06*e*e - 1.85810e-03*e + 6.16434e-01);                       
                            else t_pass3 = t_pass2 - (1.09187e-06*700*700 - 1.85810e-03*700 + 6.16434e-01);                       

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
                case 3: {
                            if(e < 700.) t_pass3 = t_pass2 - (-5.14276e-06*e*e + 4.32677e-03*e - 8.73504e-01);                        
                            else t_pass3 = t_pass2 - (-5.14276e-06*700*700 + 4.32677e-03*700 - 8.73504e-01);                        

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
                case 4: {
                            if(e < 650.) t_pass3 = t_pass2 - (1.93255e-06*e*e - 3.17103e-03*e + 7.11106e-01);                       
                            else t_pass3 = t_pass2 - (1.93255e-06*650*650 - 3.17103e-03*650 + 7.11106e-01);                       

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
                case 5: {
                            if(e < 650.) t_pass3 = t_pass2 - (-2.52727e-07*e*e - 1.66121e-03*e + 4.97491e-01);                        
                            else t_pass3 = t_pass2 - (-2.52727e-07*650*650 - 1.66121e-03*650 + 4.97491e-01);                        

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
                case 6: {
                            if(e < 700.) t_pass3 = t_pass2 - (4.57446e-06*e*e - 5.24988e-03*e + 1.41012e+00);                       
                            else t_pass3 = t_pass2 - (4.57446e-06*700*700 - 5.24988e-03*700 + 1.41012e+00);                       

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
                case 7: {
                            if(e < 700.) t_pass3 = t_pass2 - (-5.90292e-06*e*e + 5.10021e-03*e - 1.06100e+00);                        
                            else t_pass3 = t_pass2 - (-5.90292e-06*700*700 + 5.10021e-03*700 - 1.06100e+00);                        

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot[slot]->Fill(t_pass4);
                            break;
                        }
            }

            if(ph1_f1 < 0.2) t_pass4 = t_pass3 - (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
            else t_pass4 = t_pass3 - (-7.23778e-02);

            if(t_pass4 > -999.) {
                h_t_pass4->Fill(t_pass4);
                if(fabs(ph_eta1) < 1.37) h_t_pass4_emb->Fill(t_pass4);
            }
            if(t_pass3 > -999.) {
                h_t_pass3->Fill(t_pass3);
                if(fabs(ph_eta1) < 1.37) h_t_pass3_emb->Fill(t_pass3);
            }
            p_t_pass4_vs_del_eta->Fill(ph1_etas2, t_pass4);
            p_t_pass4_vs_del_phi->Fill(ph1_phis2, t_pass4);
            p_t_pass4_vs_f1->Fill(ph1_f1, t_pass4);
            p_t_pass4_vs_f3->Fill(ph1_f3, t_pass4);
            p_t_pass3_vs_del_eta->Fill(ph1_etas2, t_pass3);
            p_t_pass3_vs_del_phi->Fill(ph1_phis2, t_pass3);
            p_t_pass3_vs_f1->Fill(ph1_f1, t_pass3);
            p_t_pass3_vs_f3->Fill(ph1_f3, t_pass3);

        }
    }//END PASS 4 LOOP

    nbytes = 0;
    nb = 0;
    //Validation LOOP
    for (Long64_t jentry=0; jentry<maxentries;jentry++) {//EVENT LOOP
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( (jentry > 0) && (jentry%100000 == 0) ) cout << "Validation: Events Processed: " << jentry << endl;
        if(badOnlidCheck(ph1_maxEcell_onlId)) continue;

        if((ph1_maxEcell_gain == 2) && (m_met > 100.) && (m_met < 200.) && (ph_pt1 > 150.) && (n_ph == 1)){

            t_pass0 = tofCorrection(ph1_maxEcell_x, ph1_maxEcell_y, ph1_maxEcell_z, ph1_maxEcell_time, PV_z); 

            if(fabs(ph_eta1) < 1.37){
                if( ph_eta1 < 0){//EMBA
                    //t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emba->GetBinContent(runNumber - 297690 + 1));
                    t_pass1 = t_pass0 - ((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_emba")->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMBC
                    //t_pass1 = t_pass0 - (p_t_pass0_vs_runN_embc->GetBinContent(runNumber - 297690 + 1));
                    t_pass1 = t_pass0 - ((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_embc")->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
            }//BARREL

            if((fabs(ph_eta1) > 1.52) && (fabs(ph_eta1) < 2.5)){
                if( ph_eta1 < 0){//EMECA
                    //t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emeca->GetBinContent(runNumber - 297690 + 1));
                    t_pass1 = t_pass0 - ((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_emeca")->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
                if(ph_eta1 > 0){//EMECC
                    //t_pass1 = t_pass0 - (p_t_pass0_vs_runN_emecc->GetBinContent(runNumber - 297690 + 1));
                    t_pass1 = t_pass0 - ((TProfile*)corrections_file->Get("p_t_pass0_vs_runN_emecc")->GetBinContent(runNumber - 297690 + 1));
                    if(fabs(t_pass1) > 4.) continue;
                }
            }//ENDCAP

            onlId = ph1_maxEcell_onlId>>32; 
            febN = getFebNo(onlId);
            //t_pass2 = t_pass1 - (p_t_pass1_vs_febN->GetBinContent(febN));
            t_pass2 = t_pass1 - ((TProfile*)corrections_file->Get("p_t_pass1_vs_febN")->GetBinContent(febN));

            e = ph1_maxEcell_E/1000;
            slot = (int)(febN/32);
            t_pass3 = -999.;
            t_pass4 = -999.;
            switch(slot){
                case 0: {
                            if(e < 650.) t_pass3 = t_pass2 - (1.34131e-06*e*e - 2.59592e-03*e + 5.98932e-01);                       
                            else t_pass3 = t_pass2 - (1.34131e-06*650*650 - 2.59592e-03*650 + 5.98932e-01);                       

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
                case 1: {
                            if(e < 650.) t_pass3 = t_pass2 - (2.30630e-06*e*e - 3.59088e-03*e + 8.37710e-01);                       
                            else t_pass3 = t_pass2 - (2.30630e-06*650*650 - 3.59088e-03*650 + 8.37710e-01);                       

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
                case 2: {
                            if(e < 700.) t_pass3 = t_pass2 - (1.09187e-06*e*e - 1.85810e-03*e + 6.16434e-01);                       
                            else t_pass3 = t_pass2 - (1.09187e-06*700*700 - 1.85810e-03*700 + 6.16434e-01);                       

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
                case 3: {
                            if(e < 700.) t_pass3 = t_pass2 - (-5.14276e-06*e*e + 4.32677e-03*e - 8.73504e-01);                        
                            else t_pass3 = t_pass2 - (-5.14276e-06*700*700 + 4.32677e-03*700 - 8.73504e-01);                        

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
                case 4: {
                            if(e < 650.) t_pass3 = t_pass2 - (1.93255e-06*e*e - 3.17103e-03*e + 7.11106e-01);                       
                            else t_pass3 = t_pass2 - (1.93255e-06*650*650 - 3.17103e-03*650 + 7.11106e-01);                       

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
                case 5: {
                            if(e < 650.) t_pass3 = t_pass2 - (-2.52727e-07*e*e - 1.66121e-03*e + 4.97491e-01);                        
                            else t_pass3 = t_pass2 - (-2.52727e-07*650*650 - 1.66121e-03*650 + 4.97491e-01);                        

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
                case 6: {
                            if(e < 700.) t_pass3 = t_pass2 - (4.57446e-06*e*e - 5.24988e-03*e + 1.41012e+00);                       
                            else t_pass3 = t_pass2 - (4.57446e-06*700*700 - 5.24988e-03*700 + 1.41012e+00);                       

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
                case 7: {
                            if(e < 700.) t_pass3 = t_pass2 - (-5.90292e-06*e*e + 5.10021e-03*e - 1.06100e+00);                        
                            else t_pass3 = t_pass2 - (-5.90292e-06*700*700 + 5.10021e-03*700 - 1.06100e+00);                        

                            h_t_pass3_slot_val[slot]->Fill(t_pass3);
                            p_t_pass3_vs_E_cell_slot_val[slot]->Fill(e, t_pass3);
                            p_t_pass3_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass3);
                            p_t_pass3_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass3);

                            t_pass4 = t_pass3;

                            if(ph1_f1 < 0.2) t_pass4 -= (8.77451*ph1_f1*ph1_f1 - 2.66289*ph1_f1 + 1.71838e-01);
                            else t_pass4 -= (-7.23778e-02);

                            h_t_pass4_slot_val[slot]->Fill(t_pass4);
                            p_t_pass4_vs_channel_no_slot_val[slot]->Fill(getChannelNo(onlId)%128, t_pass4);
                            p_t_pass4_vs_channel_no_sym_slot_val[slot%4]->Fill(getChannelNo(onlId)%128, t_pass4);
                            break;
                        }
            }
            if(t_pass3 > -999.) h_t_pass3_val->Fill(t_pass3);
            if(fabs(ph_eta1) < 1.37) h_t_pass3_emb_val->Fill(t_pass3);
            if(t_pass4 > -999.) {
                h_t_pass4_val->Fill(t_pass4);
                if(fabs(ph_eta1) < 1.37) h_t_pass4_emb_val->Fill(t_pass4);
            }
        }
    }//END VALIDATION LOOP

    nbytes = 0;
    nb = 0;

    writeHists();

}

void CalibrationAnalysis::bookHists(){

    h_t_uncorr = new TH1F("h_t_uncorr", "t_{#gamma}^{uncorrected}", 10000, -200., 200.);
    h_t_pass0 = new TH1F("h_t_pass0", "t_{#gamma}^{pass0}", 10000, -200., 200.);
    h_t_pass1 = new TH1F("h_t_pass1", "t_{#gamma}^{pass1}", 10000, -200., 200.);
    h_t_pass1_vs_e = new TH2F("h_t_pass1_vs_e", "t_{#gamma}^{pass1} vs E_{#gamma}", 1000, 150., 2000., 10000, -200., 200.);
    h_t_pass1_vs_et = new TH2F("h_t_pass1_vs_et", "t_{#gamma}^{pass1} vs E_{T,#gamma}", 1000, 150., 1000., 10000, -200., 200.);
    h_t_pass0_emb = new TH1F("h_t_pass0_emb", "t_{#gamma}^{pass0}, emb", 10000, -200., 200.);
    h_t_pass1_emb = new TH1F("h_t_pass1_emb", "t_{#gamma}^{pass1}, emb", 10000, -200., 200.);
    h_t_pass0_emec = new TH1F("h_t_pass0_emec", "t_{#gamma}^{pass0}, emec", 10000, -200., 200.);
    h_t_pass1_emec = new TH1F("h_t_pass1_emec", "t_{#gamma}^{pass1}, emec", 10000, -200., 200.);
    h_t_pass0_emba = new TH1F("h_t_pass0_emba", "t_{#gamma}^{pass0}, emba", 10000, -200., 200.);
    h_t_pass1_emba = new TH1F("h_t_pass1_emba", "t_{#gamma}^{pass1}, emba", 10000, -200., 200.);
    h_t_pass0_emeca = new TH1F("h_t_pass0_emeca", "t_{#gamma}^{pass0}, emeca", 10000, -200., 200.);
    h_t_pass1_emeca = new TH1F("h_t_pass1_emeca", "t_{#gamma}^{pass1}, emeca", 10000, -200., 200.);
    h_t_pass0_embc = new TH1F("h_t_pass0_embc", "t_{#gamma}^{pass0}, embc", 10000, -200., 200.);
    h_t_pass1_embc = new TH1F("h_t_pass1_embc", "t_{#gamma}^{pass1}, embc", 10000, -200., 200.);
    h_t_pass0_emecc = new TH1F("h_t_pass0_emecc", "t_{#gamma}^{pass0}, emecc", 10000, -200., 200.);
    h_t_pass1_emecc = new TH1F("h_t_pass1_emecc", "t_{#gamma}^{pass1}, emecc", 10000, -200., 200.);
    p_t_pass0_vs_runN = new TProfile("p_t_pass0_vs_runN", "Profile of t_{#gamma}^{pass0} vs Run Number", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN = new TProfile("p_t_pass1_vs_runN", "Profile of t_{#gamma}^{pass1} vs Run Number", 13830, 297690, 311520, -10., 10.);
    p_t_pass0_vs_runN_emb = new TProfile("p_t_pass0_vs_runN_emb", "Profile of t_{#gamma}^{pass0} vs Run Number, emb", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emb = new TProfile("p_t_pass1_vs_runN_emb", "Profile of t_{#gamma}^{pass1} vs Run Number, emb", 13830, 297690, 311520, -10., 10.);
    p_t_pass0_vs_runN_emec = new TProfile("p_t_pass0_vs_runN_emec", "Profile of t_{#gamma}^{pass0} vs Run Number, emec", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emec = new TProfile("p_t_pass1_vs_runN_emec", "Profile of t_{#gamma}^{pass1} vs Run Number, emec", 13830, 297690, 311520, -10., 10.);
    p_t_pass0_vs_runN_emba = new TProfile("p_t_pass0_vs_runN_emba", "Profile of t_{#gamma}^{pass0} vs Run Number, emba", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emba = new TProfile("p_t_pass1_vs_runN_emba", "Profile of t_{#gamma}^{pass1} vs Run Number, emba", 13830, 297690, 311520, -10., 10.);
    p_t_pass0_vs_runN_emeca = new TProfile("p_t_pass0_vs_runN_emeca", "Profile of t_{#gamma}^{pass0} vs Run Number, emeca", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emeca = new TProfile("p_t_pass1_vs_runN_emeca", "Profile of t_{#gamma}^{pass1} vs Run Number, emeca", 13830, 297690, 311520, -10., 10.);
    p_t_pass0_vs_runN_embc = new TProfile("p_t_pass0_vs_runN_embc", "Profile of t_{#gamma}^{pass0} vs Run Number, embc", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_embc = new TProfile("p_t_pass1_vs_runN_embc", "Profile of t_{#gamma}^{pass1} vs Run Number, embc", 13830, 297690, 311520, -10., 10.);
    p_t_pass0_vs_runN_emecc = new TProfile("p_t_pass0_vs_runN_emecc", "Profile of t_{#gamma}^{pass0} vs Run Number, emecc", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emecc = new TProfile("p_t_pass1_vs_runN_emecc", "Profile of t_{#gamma}^{pass1} vs Run Number, emecc", 13830, 297690, 311520, -10., 10.);

    p_t_pass1_vs_e_emba = new TProfile("p_t_pass1_vs_e_emba", "Profile of t_{#gamma}^{pass1} vs Energy, emba", 100, 100, 2100, -10., 10.);
    p_t_pass1_vs_e_embc = new TProfile("p_t_pass1_vs_e_embc", "Profile of t_{#gamma}^{pass1} vs Energy, embc", 100, 100, 2100, -10., 10.);
    p_t_pass1_vs_e_emeca = new TProfile("p_t_pass1_vs_e_emeca", "Profile of t_{#gamma}^{pass1} vs Energy, emeca", 100, 100, 2100, -10., 10.);
    p_t_pass1_vs_e_emecc = new TProfile("p_t_pass1_vs_e_emecc", "Profile of t_{#gamma}^{pass1} vs Energy, emecc", 100, 100, 2100, -10., 10.);

    p_t_pass1_vs_febN = new TProfile("p_t_pass1_vs_febN", "Profile of t_{#gamma}^{pass1} vs FEB Number", 620, 1, 620, -10., 10.);
    p_t_pass2_vs_febN = new TProfile("p_t_pass2_vs_febN", "Profile of t_{#gamma}^{pass1} vs FEB Number", 620, 1, 620, -10., 10.);
    h_t_pass2 = new TH1F("h_t_pass2", "t_{#gamma}^{pass1}", 10000, -200., 200.);

    h_t_pass1_val = new TH1F("h_t_pass1_val", "t_{#gamma}^{pass1}", 10000, -200., 200.);
    p_t_pass1_vs_febN_val = new TProfile("p_t_pass1_vs_febN_val", "Profile of t_{#gamma}^{pass1} vs FEB Number", 620, 1, 620, -10., 10.);
    p_t_pass2_vs_febN_val = new TProfile("p_t_pass2_vs_febN_val", "Profile of t_{#gamma}^{pass1} vs FEB Number, Validation", 620, 1, 620, -10., 10.);
    h_t_pass2_val = new TH1F("h_t_pass2_val", "t_{#gamma}^{pass1}, Validation", 10000, -200., 200.);

    h_t_pass2_narrow = new TH1F("h_t_pass2_narrow", "t_{#gamma}^{pass1}, Validation", 10000, -200., 200.);
    h_t_pass2_wide = new TH1F("h_t_pass2_wide", "t_{#gamma}^{pass1}, Validation", 10000, -200., 200.);

    h_t_uncorr_val = new TH1F("h_t_uncorr_val", "t_{#gamma}^{uncorrected}", 10000, -200., 200.);
    h_t_pass1_val = new TH1F("h_t_pass1_val", "t_{#gamma}^{pass1}", 10000, -200., 200.);
    h_t_pass1_emb_val = new TH1F("h_t_pass1_emb_val", "t_{#gamma}^{pass1}, emb", 10000, -200., 200.);
    h_t_pass1_emec_val = new TH1F("h_t_pass1_emec_val", "t_{#gamma}^{pass1}, emec", 10000, -200., 200.);
    h_t_pass1_emba_val = new TH1F("h_t_pass1_emba_val", "t_{#gamma}^{pass1}, emba", 10000, -200., 200.);
    h_t_pass1_emeca_val = new TH1F("h_t_pass1_emeca_val", "t_{#gamma}^{pass1}, emeca", 10000, -200., 200.);
    h_t_pass1_embc_val = new TH1F("h_t_pass1_embc_val", "t_{#gamma}^{pass1}, embc", 10000, -200., 200.);
    h_t_pass1_emecc_val = new TH1F("h_t_pass1_emecc_val", "t_{#gamma}^{pass1}, emecc", 10000, -200., 200.);
    p_t_pass1_vs_runN_val = new TProfile("p_t_pass1_vs_runN_val", "Profile of t_{#gamma}^{pass1} vs Run Number", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emb_val = new TProfile("p_t_pass1_vs_runN_emb_val", "Profile of t_{#gamma}^{pass1} vs Run Number, emb", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emec_val = new TProfile("p_t_pass1_vs_runN_emec_val", "Profile of t_{#gamma}^{pass1} vs Run Number, emec", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emba_val = new TProfile("p_t_pass1_vs_runN_emba_val", "Profile of t_{#gamma}^{pass1} vs Run Number, emba", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emeca_val = new TProfile("p_t_pass1_vs_runN_emeca_val", "Profile of t_{#gamma}^{pass1} vs Run Number, emeca", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_embc_val = new TProfile("p_t_pass1_vs_runN_embc_val", "Profile of t_{#gamma}^{pass1} vs Run Number, embc", 13830, 297690, 311520, -10., 10.);
    p_t_pass1_vs_runN_emecc_val = new TProfile("p_t_pass1_vs_runN_emeca_val", "Profile of t_{#gamma}^{pass1} vs Run Number, emeca", 13830, 297690, 311520, -10., 10.);

    //Good timing
    h_t_good = new TH1F("h_t_good", "Good Timing", 10000, -100, 100);
    h_e_good = new TH1F("h_e_good", "Energy, good timing", 200, 0, 2000);
    h_et_good = new TH1F("h_et_good", "Transverse Energy, good timing", 200, 0, 2000);
    h_met_good = new TH1F("h_met_good", "Missing Transverse Energy, good timing", 200, 0, 2000);
    h_Ecell_E_good = new TH1F("h_Ecell_E_good", "Max Ecell E, good timing", 200, 0, 2000);
    h_Ecell_x_good = new TH1F("h_Ecell_x_good", "Max Ecell x, good timing", 200, -1000, 1000);
    h_Ecell_y_good = new TH1F("h_Ecell_y_good", "May Ecell y, good timing", 200, -1000, 1000);
    h_Ecell_z_good = new TH1F("h_Ecell_z_good", "Maz Ecell z, good timing", 200, -1000, 1000);
    h_eta_good = new TH1F("h_eta_good", "Eta, good timing", 200, -2.5, 2.5);
    h_phi_good = new TH1F("h_phi_good", "Phi, good timing", 200, -3.15, 3.15);
    h_cl_eta_good = new TH1F("h_cl_eta_good", "Cluster Eta, good timing", 200, -2.5, 2.5);
    h_cl_phi_good = new TH1F("h_cl_phi_good", "Cluster Phi, good timing", 200, -3.15, 3.15);
    h_f1_good = new TH1F("h_f1_good", "f1, good timing", 200, 0, 1);
    h_f3_good = new TH1F("h_f3_good", "f3, good timing", 200, 0, 1);

    //Bad timing
    h_t_bad = new TH1F("h_t_bad", "Good Timing", 10000, -100, 100);
    h_e_bad = new TH1F("h_e_bad", "Energy, bad timing", 200, 0, 2000);
    h_et_bad = new TH1F("h_et_bad", "Transverse Energy, bad timing", 200, 0, 2000);
    h_met_bad = new TH1F("h_met_bad", "Missing Transverse Energy, bad timing", 200, 0, 2000);
    h_Ecell_E_bad = new TH1F("h_Ecell_E_bad", "Max Ecell E, bad timing", 200, 0, 2000);
    h_Ecell_x_bad = new TH1F("h_Ecell_x_bad", "Max Ecell x, bad timing", 200, -1000, 1000);
    h_Ecell_y_bad = new TH1F("h_Ecell_y_bad", "May Ecell y, bad timing", 200, -1000, 1000);
    h_Ecell_z_bad = new TH1F("h_Ecell_z_bad", "Maz Ecell z, bad timing", 200, -1000, 1000);
    h_eta_bad = new TH1F("h_eta_bad", "Eta, bad timing", 200, -2.5, 2.5);
    h_phi_bad = new TH1F("h_phi_bad", "Phi, bad timing", 200, -3.15, 3.15);
    h_cl_eta_bad = new TH1F("h_cl_eta_bad", "Cluster Eta, bad timing", 200, -2.5, 2.5);
    h_cl_phi_bad = new TH1F("h_cl_phi_bad", "Cluster Phi, bad timing", 200, -3.15, 3.15);
    h_f1_bad = new TH1F("h_f1_bad", "f1, bad timing", 200, 0, 1);
    h_f3_bad = new TH1F("h_f3_bad", "f3, bad timing", 200, 0, 1);

    h_t_bad_stdm2 = new TH1F("h_t_bad_stdm2", "Events from bad stdm2 timing", 10000, -100, 100);

    h_t_32febs = new TH1F("h_t_32febs", "t_{#gamma}^{pass1}, First 32 febs", 100, -10., 10.);
    h_t_corr_32febs = new TH1F("h_t_corr_32febs", "t_{#gamma}^{pass1}, corrected, First 32 febs", 100, -10., 10.);

    h_t_32febs_wide = new TH1F("h_t_32febs_wide", "t_{#gamma}^{pass1}, Febs 131-162", 100, -10., 10.);
    h_t_corr_32febs_wide = new TH1F("h_t_corr_32febs_wide", "t_{#gamma}^{pass1}, corrected, Febs 131-162", 100, -10., 10.);

    for(int i = 1; i < 257; i++){
        h_t_fn[i-1] = new TH1F(Form("h_t_fn_%d", i), Form("t_{#gamma}^{pass1}, fn_%d", i), 100, -10., 10.);
    }

    for(int i = 0; i < 8; i++){
        h_t_febs_slot[i] = new TH1F(Form("h_t_febs_slot_%d", i), Form("t_{#gamma}^{pass1}, slot %d", i), 100, -10., 10.);
        h_t_corr_febs_slot[i] = new TH1F(Form("h_t_corr_febs_slot_%d", i), Form("t_{#gamma}^{pass1}, slot %d", i), 100, -10., 10.);
        h_e_slot[i] = new TH1F(Form("h_e_slot_%d", i), Form("Energy, slot %d", i), 200, 0, 2000);
        h_e_cell_slot[i] = new TH1F(Form("h_e_cell_slot_%d", i), Form("Cell Energy, slot %d", i), 200, 0, 2000);
        h_et_slot[i] = new TH1F(Form("h_et_slot_%d", i), Form("Transverse Energy, slot %d", i), 200, 0, 2000);
        h_eta_slot[i] = new TH1F(Form("h_eta_slot_%d", i), Form("Eta, slot %d", i), 200, -2.5, 2.5);
        h_phi_slot[i] = new TH1F(Form("h_phi_slot_%d", i), Form("Phi, slot %d", i), 200, -3.15, 3.15);
        p_t_vs_E_slot[i] = new TProfile(Form("p_t_vs_E_slot_%d", i), Form("Profile of t vs E, slot_%d", i), 100, 100, 2100, -10., 10.);
        p_t_corr_vs_E_slot[i] = new TProfile(Form("p_t_corr_vs_E_slot_%d", i), Form("Profile of t corrected vs E, slot_%d", i), 100, 100, 2100, -10., 10.);
        p_t_vs_ET_slot[i] = new TProfile(Form("p_t_vs_ET_slot_%d", i), Form("Profile of t vs ET, slot_%d", i), 100, 100, 2100, -10., 10.);
        p_t_corr_vs_ET_slot[i] = new TProfile(Form("p_t_corr_vs_ET_slot_%d", i), Form("Profile of t corrected vs ET, slot_%d", i), 100, 100, 2100, -10., 10.);
        p_t_vs_E_cell_slot[i] = new TProfile(Form("p_t_vs_E_cell_slot_%d", i), Form("Profile of t vs E_cell, slot_%d", i), 100, 100, 2100, -10., 10.);
        p_t_pass2_vs_E_cell_slot[i] = new TProfile(Form("p_t_pass2_vs_E_cell_slot_%d", i), Form("Profile of t corrected vs E_cell, slot_%d", i), 100, 100, 2100, -10., 10.);
        h_t_pass3_slot[i] = new TH1F(Form("h_t_pass3_slot_%d", i), Form("t_{#gamma}^{pass3}, slot %d", i), 100, -10., 10.);
        p_t_pass3_vs_E_cell_slot[i] = new TProfile(Form("p_t_pass3_vs_E_cell_slot_%d", i), Form("Profile of t_{#gamma}^{pass3} vs E_cell, slot_%d", i), 100, 100, 2100, -10., 10.);
        p_t_pass3_vs_channel_no_slot[i] = new TProfile(Form("p_t_pass3_vs_channel_no_slot_%d", i), Form("Profile of t_{#gamma}^{pass3} vs Channel Number, slot_%d", i), 128, 0, 127, -10., 10.);
        h_t_pass4_slot[i] = new TH1F(Form("h_t_pass4_slot_%d", i), Form("t_{#gamma}^{pass4}, slot %d", i), 100, -10., 10.);
        p_t_pass4_vs_channel_no_slot[i] = new TProfile(Form("p_t_pass4_vs_channel_no_slot_%d", i), Form("Profile of t_{#gamma}^{pass4} vs Channel Number, slot_%d", i), 128, 0, 127, -10., 10.);
        h_t_pass3_slot_val[i] = new TH1F(Form("h_t_pass3_slot_val_%d", i), Form("t_{#gamma}^{pass3}, slot_val %d", i), 100, -10., 10.);
        p_t_pass3_vs_E_cell_slot_val[i] = new TProfile(Form("p_t_pass3_vs_E_cell_slot_val_%d", i), Form("Profile of t_{#gamma}^{pass3} vs E_cell, slot_val_%d", i), 100, 100, 2100, -10., 10.);
        p_t_pass3_vs_channel_no_slot_val[i] = new TProfile(Form("p_t_pass3_vs_channel_no_slot_val_%d", i), Form("Profile of t_{#gamma}^{pass3} vs Channel Number, slot_val_%d", i), 128, 0, 127, -10., 10.);
        h_t_pass4_slot_val[i] = new TH1F(Form("h_t_pass4_slot_val_%d", i), Form("t_{#gamma}^{pass4}, slot_val %d", i), 100, -10., 10.);
        p_t_pass4_vs_channel_no_slot_val[i] = new TProfile(Form("p_t_pass4_vs_channel_no_slot_val_%d", i), Form("Profile of t_{#gamma}^{pass4} vs Channel Number, slot_val_%d", i), 128, 0, 127, -10., 10.);

    }

    for(int i = 0; i < 4; i++){
        p_t_pass3_vs_channel_no_sym_slot[i] = new TProfile(Form("p_t_pass3_vs_channel_no_sym_slot_%d", i), Form("Profile of t_{#gamma}^{pass3} vs Channel Number, sym_slot_%d", i), 128, 0, 127, -10., 10.);
        p_t_pass4_vs_channel_no_sym_slot[i] = new TProfile(Form("p_t_pass4_vs_channel_no_sym_slot_%d", i), Form("Profile of t_{#gamma}^{pass4} vs Channel Number, sym_slot_%d", i), 128, 0, 127, -10., 10.);
        p_t_pass3_vs_channel_no_sym_slot_val[i] = new TProfile(Form("p_t_pass3_vs_channel_no_sym_slot_val_%d", i), Form("Profile of t_{#gamma}^{pass3} vs Channel Number, sym_slot_val_%d", i), 128, 0, 127, -10., 10.);
        p_t_pass4_vs_channel_no_sym_slot_val[i] = new TProfile(Form("p_t_pass4_vs_channel_no_sym_slot_val_%d", i), Form("Profile of t_{#gamma}^{pass4} vs Channel Number, sym_slot_val_%d", i), 128, 0, 127, -10., 10.);
    }

    h_t_pass3 = new TH1F("h_t_pass3", "t_{#gamma}^{pass3}", 100, -10., 10.);
    h_t_pass3_emb = new TH1F("h_t_pass3_emb", "t_{#gamma}^{pass3}, emb", 100, -10., 10.);
    h_t_pass4 = new TH1F("h_t_pass4", "t_{#gamma}^{pass4}", 100, -10., 10.);
    h_t_pass4_emb = new TH1F("h_t_pass4_emb", "t_{#gamma}^{pass4}, emb", 100, -10., 10.);

    h_t_pass3_val = new TH1F("h_t_pass3_val", "t_{#gamma}^{pass3}_val", 100, -10., 10.);
    h_t_pass3_emb_val = new TH1F("h_t_pass3_emb_val", "t_{#gamma}^{pass3}, emb_val", 100, -10., 10.);
    h_t_pass4_val = new TH1F("h_t_pass4_val", "t_{#gamma}^{pass4}_val", 100, -10., 10.);
    h_t_pass4_emb_val = new TH1F("h_t_pass4_emb_val", "t_{#gamma}^{pass4}, emb_val", 100, -10., 10.);

    p_t_pass3_vs_del_eta = new TProfile("p_t_pass3_vs_del_eta", "Profile of t_{#gamma}^{pass3} vs #delta#eta", 100, -0.6, 0.6, -10., 10.);
    p_t_pass3_vs_del_phi = new TProfile("p_t_pass3_vs_del_phi", "Profile of t_{#gamma}^{pass3} vs #delta#phi", 100, -1., 1., -10., 10.);
    p_t_pass3_vs_f1 = new TProfile("p_t_pass3_vs_f1", "Profile of t_{#gamma}^{pass3} vs f1", 100, 0, 0.6, -10., 10.);
    p_t_pass3_vs_f3 = new TProfile("p_t_pass3_vs_f3", "Profile of t_{#gamma}^{pass3} vs f3", 100, 0, 0.15, -10., 10.);

    p_t_pass4_vs_del_eta = new TProfile("p_t_pass4_vs_del_eta", "Profile of t_{#gamma}^{pass4} vs #delta#eta", 100, -0.6, 0.6, -10., 10.);
    p_t_pass4_vs_del_phi = new TProfile("p_t_pass4_vs_del_phi", "Profile of t_{#gamma}^{pass4} vs #delta#phi", 100, -1., 1., -10., 10.);
    p_t_pass4_vs_f1 = new TProfile("p_t_pass4_vs_f1", "Profile of t_{#gamma}^{pass4} vs f1", 100, 0, 0.6, -10., 10.);
    p_t_pass4_vs_f3 = new TProfile("p_t_pass4_vs_f3", "Profile of t_{#gamma}^{pass4} vs f3", 100, 0, 0.15, -10., 10.);

    h_t_febs_narrow = new TH1F("h_t_febs_narrow", "t_{#gamma}^{pass1}, narrow FEBs", 100, -10., 10.);
    h_t_corr_febs_narrow = new TH1F("h_t_corr_febs_narrow", "t_{#gamma}^{pass1}, narrow FEBs", 100, -10., 10.);
    h_e_narrow = new TH1F("h_e_narrow", "Energy, narrow FEBs", 200, 0, 2000);
    h_e_cell_narrow = new TH1F("h_e_cell_narrow", "Cell Energy, narrow FEBs", 200, 0, 2000);
    h_et_narrow = new TH1F("h_et_narrow", "Transverse Energy, narrow FEBs", 200, 0, 2000);
    h_eta_narrow = new TH1F("h_eta_narrow", "Eta, narrow FEBs", 200, -2.5, 2.5);
    h_phi_narrow = new TH1F("h_phi_narrow", "Phi, narrow FEBs", 200, -3.15, 3.15);

    h_t_febs_wide = new TH1F("h_t_febs_wide", "t_{#gamma}^{pass1}, wide FEBs", 100, -10., 10.);
    h_t_corr_febs_wide = new TH1F("h_t_corr_febs_wide", "t_{#gamma}^{pass1}, wide FEBs", 100, -10., 10.);
    h_e_wide = new TH1F("h_e_wide", "Energy, wide FEBs", 200, 0, 2000);
    h_e_cell_wide = new TH1F("h_e_cell_wide", "Cell Energy, wide FEBs", 200, 0, 2000);
    h_et_wide = new TH1F("h_et_wide", "Transverse Energy, wide FEBs", 200, 0, 2000);
    h_eta_wide = new TH1F("h_eta_wide", "Eta, wide FEBs", 200, -2.5, 2.5);
    h_phi_wide = new TH1F("h_phi_wide", "Phi, wide FEBs", 200, -3.15, 3.15);

    h_t_febs_emec = new TH1F("h_t_febs_emec", "t_{#gamma}^{pass1}, emec FEBs", 100, -10., 10.);
    h_t_corr_febs_emec = new TH1F("h_t_corr_febs_emec", "t_{#gamma}^{pass1}, emec FEBs", 100, -10., 10.);
    h_e_emec = new TH1F("h_e_emec", "Energy, emec FEBs", 200, 0, 2000);
    h_e_cell_emec = new TH1F("h_e_cell_emec", "Cell Energy, emec FEBs", 200, 0, 2000);
    h_et_emec = new TH1F("h_et_emec", "Transverse Energy, emec FEBs", 200, 0, 2000);
    h_eta_emec = new TH1F("h_eta_emec", "Eta, emec FEBs", 200, -2.5, 2.5);
    h_phi_emec = new TH1F("h_phi_emec", "Phi, emec FEBs", 200, -3.15, 3.15);
    p_t_corr_emec_vs_subslot = new TProfile("p_t_corr_emec_vs_subslot", "Profile of t_{#gamma}^{pass1} vs Subslot Number, EMEC", 7, 22, 28, -10., 10.);

    h_e_fullFEB = new TH1F("h_e_fullFEB", "Energy", 200, 0, 2000);

}

void CalibrationAnalysis::writeHists(){

    TFile f("NtuplePlots/calibration_exot6.root", "recreate");

    p_t_pass0_vs_runN_emba->Write();
    p_t_pass0_vs_runN_emeca->Write();
    p_t_pass0_vs_runN_embc->Write();
    p_t_pass0_vs_runN_emecc->Write();
    p_t_pass1_vs_febN->Write();

    h_t_uncorr->Write();
    h_t_pass0->Write();
    h_t_pass1->Write();
    h_t_pass1_vs_e->Write();
    h_t_pass1_vs_et->Write();
    p_t_pass0_vs_runN->Write();
    p_t_pass1_vs_runN->Write();
    h_t_pass0_emb->Write();
    h_t_pass1_emb->Write();
    h_t_pass0_emec->Write();
    h_t_pass1_emec->Write();
    p_t_pass0_vs_runN_emb->Write();
    p_t_pass1_vs_runN_emb->Write();
    p_t_pass0_vs_runN_emec->Write();
    p_t_pass1_vs_runN_emec->Write();
    h_t_pass0_emba->Write();
    h_t_pass1_emba->Write();
    h_t_pass0_emeca->Write();
    h_t_pass1_emeca->Write();
    p_t_pass1_vs_runN_emba->Write();
    p_t_pass1_vs_runN_emeca->Write();
    h_t_pass0_embc->Write();
    h_t_pass1_embc->Write();
    h_t_pass0_emecc->Write();
    h_t_pass1_emecc->Write();
    p_t_pass1_vs_runN_embc->Write();
    p_t_pass1_vs_runN_emecc->Write();

    p_t_pass1_vs_e_emba->Write();
    p_t_pass1_vs_e_embc->Write();
    p_t_pass1_vs_e_emeca->Write();
    p_t_pass1_vs_e_emecc->Write();

    p_t_pass2_vs_febN->Write();
    h_t_pass2->Write();

    h_t_pass1_val->Write();
    p_t_pass1_vs_febN_val->Write();
    p_t_pass2_vs_febN_val->Write();
    h_t_pass2_val->Write();

    h_t_pass2_narrow->Write();
    h_t_pass2_wide->Write();

    h_t_uncorr_val->Write();
    h_t_pass1_val->Write();
    p_t_pass1_vs_runN_val->Write();
    h_t_pass1_emb_val->Write();
    h_t_pass1_emec_val->Write();
    p_t_pass1_vs_runN_emb_val->Write();
    p_t_pass1_vs_runN_emec_val->Write();
    h_t_pass1_emba_val->Write();
    h_t_pass1_emeca_val->Write();
    p_t_pass1_vs_runN_emba_val->Write();
    p_t_pass1_vs_runN_emeca_val->Write();
    h_t_pass1_embc_val->Write();
    h_t_pass1_emecc_val->Write();
    p_t_pass1_vs_runN_embc_val->Write();
    p_t_pass1_vs_runN_emecc_val->Write();

    h_t_good->Write();
    h_e_good->Write();
    h_et_good->Write();
    h_met_good->Write();
    h_Ecell_E_good->Write();
    h_Ecell_x_good->Write();
    h_Ecell_y_good->Write();
    h_Ecell_z_good->Write();
    h_eta_good->Write();
    h_phi_good->Write();
    h_cl_eta_good->Write();
    h_cl_phi_good->Write();
    h_f1_good->Write();
    h_f3_good->Write();
    h_t_bad->Write();
    h_e_bad->Write();
    h_et_bad->Write();
    h_met_bad->Write();
    h_Ecell_E_bad->Write();
    h_Ecell_x_bad->Write();
    h_Ecell_y_bad->Write();
    h_Ecell_z_bad->Write();
    h_eta_bad->Write();
    h_phi_bad->Write();
    h_cl_eta_bad->Write();
    h_cl_phi_bad->Write();
    h_f1_bad->Write();
    h_f3_bad->Write();

    h_t_bad_stdm2->Write();

    h_t_32febs->Write();
    h_t_corr_32febs->Write();

    h_t_32febs_wide->Write();
    h_t_corr_32febs_wide->Write();

    for(int i = 1; i < 257; i++){
        h_t_fn[i-1]->Write();
    }

    for(int i = 0; i < 8; i++){
        h_t_febs_slot[i]->Write();
        h_t_corr_febs_slot[i]->Write();
        h_e_slot[i]->Write();
        h_e_cell_slot[i]->Write();
        h_et_slot[i]->Write();
        h_eta_slot[i]->Write();
        h_phi_slot[i]->Write();
        p_t_vs_E_slot[i]->Write();
        p_t_corr_vs_E_slot[i]->Write();
        p_t_vs_ET_slot[i]->Write();
        p_t_corr_vs_ET_slot[i]->Write();
        p_t_vs_E_cell_slot[i]->Write();
        p_t_pass2_vs_E_cell_slot[i]->Write();
        h_t_pass3_slot[i]->Write();
        p_t_pass3_vs_E_cell_slot[i]->Write();
        p_t_pass3_vs_channel_no_slot[i]->Write();
        h_t_pass4_slot[i]->Write();
        p_t_pass4_vs_channel_no_slot[i]->Write();
        h_t_pass3_slot_val[i]->Write();
        p_t_pass3_vs_E_cell_slot_val[i]->Write();
        p_t_pass3_vs_channel_no_slot_val[i]->Write();
        h_t_pass4_slot_val[i]->Write();
        p_t_pass4_vs_channel_no_slot_val[i]->Write();
    }

    for(int i = 0; i < 4; i++){
        p_t_pass3_vs_channel_no_sym_slot[i]->Write();
        p_t_pass4_vs_channel_no_sym_slot[i]->Write();
        p_t_pass3_vs_channel_no_sym_slot_val[i]->Write();
        p_t_pass4_vs_channel_no_sym_slot_val[i]->Write();
    }

    h_t_pass3->Write();
    h_t_pass3_emb->Write();
    h_t_pass4->Write();
    h_t_pass4_emb->Write();

    h_t_pass3_val->Write();
    h_t_pass3_emb_val->Write();
    h_t_pass4_val->Write();
    h_t_pass4_emb_val->Write();

    p_t_pass3_vs_del_eta->Write();
    p_t_pass3_vs_del_phi->Write();
    p_t_pass3_vs_f1->Write();
    p_t_pass3_vs_f3->Write();

    p_t_pass4_vs_del_eta->Write();
    p_t_pass4_vs_del_phi->Write();
    p_t_pass4_vs_f1->Write();
    p_t_pass4_vs_f3->Write();

    h_t_febs_narrow->Write();
    h_t_corr_febs_narrow->Write();
    h_e_narrow->Write();
    h_e_cell_narrow->Write();
    h_et_narrow->Write();
    h_eta_narrow->Write();
    h_phi_narrow->Write();

    h_t_febs_wide->Write();
    h_t_corr_febs_wide->Write();
    h_e_wide->Write();
    h_e_cell_wide->Write();
    h_et_wide->Write();
    h_eta_wide->Write();
    h_phi_wide->Write();

    h_t_febs_emec->Write();
    h_t_corr_febs_emec->Write();
    h_e_emec->Write();
    h_e_cell_emec->Write();
    h_et_emec->Write();
    h_eta_emec->Write();
    h_phi_emec->Write();
    p_t_corr_emec_vs_subslot->Write();

    h_e_fullFEB->Write();

}


bool CalibrationAnalysis::badOnlidCheck(unsigned int onl_id){

    const unsigned int badOnlidList[2] = {950381824, 980756480}; //!

    bool badOnlid = false;
    for(int i = 0; i < 2; i++){
        if(badOnlidList[i] == onl_id) badOnlid = true;
    }
    return badOnlid;
}

int CalibrationAnalysis::getFebNo(unsigned long int onlId){
    // 620 Febs in total
    // EMB:  Sl 10-14 (5) x A/C (2) x 32FT = 320
    // EMEC: Sl 10-15 (6) x A/C (2) X 25FT = 300
    if( onlId >> 26 != 14){
        return -1;
    }

    int pn   = ( onlId >> 24 )&0x1;
    int ft   = ( onlId >> 19 )&0x1F;
    int sl   = ((onlId >> 15 )&0xF)+1;
    int benc = ( onlId >> 25 )&0x1;

    int bin = (sl - 11)*32 + ft;

    if(!pn) bin = bin + 128;

    if( sl == 10){
        bin = 256 + ft;
        if(!pn) bin = bin + 32;
    }

    if(benc){
        bin = 256 + 64 + (sl - 10)*25 + ft; // was -1
        if(!pn) bin = bin + 6*25; // was 375
    }

    return bin;
}

int CalibrationAnalysis::getSlotNo(unsigned long int onlId){

    if( onlId >> 26 != 14){
        return -1;
    }

    int pn   = ( onlId >> 24 )&0x1;
    int sl   = ((onlId >> 15 )&0xF)+1;
    int benc = ( onlId >> 25 )&0x1;

    int bin = sl - 11;

    if(!pn) bin = bin + 4;

    if( sl == 10 ){
        bin = 8;
        if(!pn) bin = 9;
    }

    if(benc){
        if( sl < 10 || sl > 15){
        }
        bin = 10 + sl -10; // was -1
        if(!pn) bin = bin + 6; // was +15
    }

    return bin;

}

int CalibrationAnalysis::getSubslotNo(unsigned long int onlId){

    int slot  = getSlotNo(onlId);   
    switch(slot){
        case 12: slot = 22;break;
        case 13: slot = 23;break;
        case 14: slot = 24;break;
        case 15: slot = 25;break;
        case 19: slot = 26;break;
        case 20: slot = 27;break;
        case 21: slot = 28;
    }
    return slot;
}

int CalibrationAnalysis::getChannelNo(unsigned long int onlId){

    // 620 Febs x 128 ch/feb = 79360 ch
    int fn = getFebNo(onlId);
    int ch = (onlId >> 8)&0x7F;
    int bin = fn * 128 + ch;  

    return bin;
}

double CalibrationAnalysis::tofCorrection(double x, double y, double z, double t_uncorr, double pv_z){

    const double c = 299.792458;

    double r_cell = sqrt( x*x + y*y + z*z);
    double r_path = sqrt( x*x + y*y + (z - pv_z)*(z - pv_z));

    double t_diff = (r_path - r_cell)/c;

    return (t_uncorr - t_diff);
}
