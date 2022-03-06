#include<string.h>
#include<stdio.h>
#include<stdlib.h>


void fit(){
//		gStyle->SetOptStat("0");	
		gStyle->SetOptFit();

		TCanvas* tmpCan = new TCanvas();
		TFile* tmpFile;
		TH1F* tmpHist;
		TH1F* metHist;
		TProfile* tmpProf;

		std::vector<string> file_names;

		file_names.push_back("../NtuplePlots/calibration_exot6.root");
		tmpFile = TFile::Open(file_names.at(0).c_str());

		tmpCan->SetLogy();

		tmpProf = (TProfile*)tmpFile->Get("p_t_pass2_vs_f1");
		TF1 *myfit = new TF1("myfit", "[0]*x*x + [1]*x + [2]", 0., 0.2);
		tmpProf->Fit("myfit", "R");
		TF1 *fit1 = tmpProf->GetFunction("myfit");	
		TF1 *myfit2 = new TF1("myfit2", "[0]", 0.2, 0.4);
		tmpProf->Fit("myfit2", "R+");
		TF1 *fit2 = tmpProf->GetFunction("myfit2");	
		tmpProf->SetTitle("Fit Result: Profile of t_{#gamma}^{pass2} vs f1");
		tmpProf->GetXaxis()->SetTitle("f1");
		tmpProf->GetYaxis()->SetTitle("Mean t_{#gamma}^{pass2} [ns]");
		tmpProf->GetYaxis()->SetRangeUser(-0.75, 0.5);
		tmpProf->Draw();
		tmpCan->SaveAs("Plots/FitResult_f1.png");

}
