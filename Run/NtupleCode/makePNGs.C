#include<string.h>
#include<stdio.h>
#include<stdlib.h>


void makePNGs(){
		//bool do1DPlots = true;
		bool do1DPlots = true;
		if(!do1DPlots) gStyle->SetOptStat("0");	
		//gStyle->SetOptFit();

		TCanvas* tmpCan = new TCanvas();
		//tmpCan->SetLeftMargin(0.15);
		//tmpCan->SetBottomMargin(0.15);
		TFile* tmpFile;

		TH1F* metHist;
		TProfile* tmpProf;
		TH1F* tmpHist;
		TStyle* style;

		std::vector<string> file_names;
		std::vector<string> iso_names;
		std::vector<string> selection_names;
		std::vector<string> region_names;
		std::vector<string> fit_names;

		std::vector<string> iso_titles;
		std::vector<string> selection_titles;
		std::vector<string> region_titles;
		std::vector<string> fit_titles;




		file_names.push_back("../NtuplePlots/calibration_exot6.root");

		//			 file_names.push_back("../NtuplePlots/1000_900_100ns_plots.root");
		//			 file_names.push_back("../NtuplePlots/1300_100_6ns_plots.root");
		//			 file_names.push_back("../NtuplePlots/1300_1200_6ns_plots.root");
		//			 file_names.push_back("../NtuplePlots/1800_100_2ns_plots.root");
		//			 file_names.push_back("../NtuplePlots/1800_1700_2ns_plots.root");
		//			 file_names.push_back("../NtuplePlots/1800_100_6ns_plots.root");
		//			 file_names.push_back("../NtuplePlots/1800_1700_6ns_plots.root");


		tmpFile = TFile::Open(file_names.at(0).c_str());

		tmpProf = (TProfile*)tmpFile->Get("p_t_pass1_corr_vs_febN");
		tmpProf->Draw();
		tmpCan->SaveAs("Plots/p_t_pass1_corr_vs_febN.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_t_pass1_corr");
		tmpHist->SetAxisRange(-10., 10.);
//		tmpCan->SetLogy();
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_pass1_corr.png");
		tmpCan->Clear();




		tmpHist = (TH1F*)tmpFile->Get("h_e_narrow");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_e_narrow.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_e_cell_narrow");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_e_cell_narrow.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_et_narrow");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_et_narrow.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_eta_narrow");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_eta_narrow.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_phi_narrow");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_phi_narrow.png");
		tmpCan->Clear();



		tmpHist = (TH1F*)tmpFile->Get("h_e_wide");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_e_wide.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_e_cell_wide");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_e_cell_wide.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_et_wide");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_et_wide.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_eta_wide");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_eta_wide.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_phi_wide");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_phi_wide.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_e_emec");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_e_emec.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_e_cell_emec");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_e_cell_emec.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_et_emec");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_et_emec.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_eta_emec");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_eta_emec.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_phi_emec");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_phi_emec.png");
		tmpCan->Clear();


		for(int i = 0; i < 8; i++){
				tmpHist = (TH1F*)tmpFile->Get(Form("h_e_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_e_slot_%d.png", i));
				tmpCan->Clear();

				tmpHist = (TH1F*)tmpFile->Get(Form("h_e_cell_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_e_cell_slot_%d.png", i));
				tmpCan->Clear();

				tmpHist = (TH1F*)tmpFile->Get(Form("h_et_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_et_slot_%d.png", i));
				tmpCan->Clear();

				tmpHist = (TH1F*)tmpFile->Get(Form("h_eta_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_eta_slot_%d.png", i));
				tmpCan->Clear();

				tmpHist = (TH1F*)tmpFile->Get(Form("h_phi_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_phi_slot_%d.png", i));
				tmpCan->Clear();


				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_vs_E_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_vs_E_slot_%d.png", i));
				tmpCan->Clear();

				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_corr_vs_E_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_corr_vs_E_slot_%d.png", i));
				tmpCan->Clear();


				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_vs_ET_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_vs_ET_slot_%d.png", i));
				tmpCan->Clear();

				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_corr_vs_ET_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_corr_vs_ET_slot_%d.png", i));
				tmpCan->Clear();


				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_vs_E_cell_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_vs_E_cell_slot_%d.png", i));
				tmpCan->Clear();

				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_corr_vs_E_cell_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_corr_vs_E_cell_slot_%d.png", i));
				tmpCan->Clear();

				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_pass2_vs_E_cell_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_pass2_vs_E_cell_slot_%d.png", i));
				tmpCan->Clear();

				tmpProf = (TProfile*)tmpFile->Get(Form("p_t_pass2_vs_channel_no_slot_%d", i));
				tmpProf->Draw();
				tmpCan->SaveAs(Form("Plots/p_t_pass2_vs_channel_no_slot_%d.png", i));
				tmpCan->Clear();


		}

		tmpCan->SetLogy();

		for(int i = 1; i < 257; i++){
				tmpHist = (TH1F*)tmpFile->Get(Form("h_t_fn_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_t_fn_%d.png", i));
				tmpCan->Clear();
		}


		tmpHist = (TH1F*)tmpFile->Get("h_t_febs_narrow");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_febs_narrow.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_t_corr_febs_narrow");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_corr_febs_narrow.png");
		tmpCan->Clear();


		tmpHist = (TH1F*)tmpFile->Get("h_t_febs_wide");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_febs_wide.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_t_corr_febs_wide");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_corr_febs_wide.png");
		tmpCan->Clear();



		tmpHist = (TH1F*)tmpFile->Get("h_t_febs_emec");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_febs_emec.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_t_corr_febs_emec");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_corr_febs_emec.png");
		tmpCan->Clear();

		tmpHist = (TH1F*)tmpFile->Get("h_t_pass2");
		tmpHist->Draw();
		tmpCan->SaveAs("Plots/h_t_pass2.png");
		tmpCan->Clear();

		for(int i = 0; i < 8; i++){
				tmpHist = (TH1F*)tmpFile->Get(Form("h_t_febs_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_t_febs_slot_%d.png", i));
				tmpCan->Clear();

				tmpHist = (TH1F*)tmpFile->Get(Form("h_t_corr_febs_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_t_corr_febs_slot_%d.png", i));
				tmpCan->Clear();

				tmpHist = (TH1F*)tmpFile->Get(Form("h_t_pass2_slot_%d", i));
				tmpHist->Draw();
				tmpCan->SaveAs(Form("Plots/h_t_pass2_slot_%d.png", i));
				tmpCan->Clear();
		}


}
