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
		if (fChain == 0) return;




		//--------------------------
		//Initialize histgrams here
		//-------------------------

		TH1F* h_MET = new TH1F("h_MET", "E_{T}^{Miss} [GeV]", 1000, 0., 1000.);







		//---------------------------------
		//Code to start the event loop here
		//---------------------------------

		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {//Actual event loop starts here
				Long64_t ientry = LoadTree(jentry);
				if (ientry < 0) break;
				nb = fChain->GetEntry(jentry);   nbytes += nb;
				// if (Cut(ientry) < 0) continue;
				if (jentry % 100000 == 0) cout << "Processed " << jentry << " events" << endl;

				//Fill Plots in this loop, like so
				h_MET->Fill(m_met/1000.);



		}//END LOOP


		//Output file
		TFile f(Form("NtuplePlots/%s_plots.root", fileName.c_str()),"recreate");

		//Write histograms to file
		h_MET->Write();


		//Close file
		f.Close();
}
