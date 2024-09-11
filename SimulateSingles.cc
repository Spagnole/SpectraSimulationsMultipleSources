#include "ReadInputData.cc"

//this function fills the simulated peak intensities
//this function requires the output from the ReadDecayScheme(), GetEfficiency() and GetPeakWidth()
//scaling paramter is used to fix for using relative efficiency and/or relative intensities
//for absolute efficiency and absolute intensities just use the default value
void FillSimulation(int source_number, double scale_int = 1.0){

		if( gammaList[source_number].size() == 0 ){
				cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
				return;
		}

		const int Nbins = hBkgr->GetXaxis()->GetNbins();
		double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
		double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
		hSimPeaks[source_number] = new TH1D(Form("hPeaks_%s", source_name.at(source_number).c_str()),Form("Simulated Source Peaks: %s",  source_name.at(source_number).c_str() ),Nbins, x_low, x_max);
		double counts;
		double width;
		for(auto MyGamma : gammaList[source_number] ){
				width = fWidth->Eval( MyGamma.gammaEn );
				counts = MyGamma.gammaInt * gEff->Eval( MyGamma.gammaEn ) * scale_int;
				for(int j = 0; j < counts; j++){
						hSimPeaks[source_number]->Fill( gRandom->Gaus(MyGamma.gammaEn,width) );
				}
		}
	hSimPeaks[source_number]->SetLineColor(kRed);
}

//this function fills the simulated single-escape peak intensities
//this function requires the output from the ReadEscapePeaks() function
//!!!USE SAME SCALING THAT WAS USED BY FillSimulation()
void FillEscapePeaks(int source_number, double scale_int = 1.0){

		const int Nbins = hBkgr->GetXaxis()->GetNbins();
		double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
		double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
		hEscPeaks[source_number] = new TH1D(Form("hEscPeaks_%s", source_name.at(source_number).c_str()),Form("Simulated Source Single-Escape Peaks: %s",  source_name.at(source_number).c_str() ),Nbins, x_low, x_max);
		hEscPeaks[source_number]->SetLineColor(kGreen+2);

		double width;
		double counts;
		for(auto MyGamma : gammaList[source_number] ){
				if( MyGamma.gammaEn < 1500. ) continue;
				counts = MyGamma.gammaInt * gEff->Eval( MyGamma.gammaEn ) * fEscPeak->Eval( MyGamma.gammaEn ) * scale_int;
				width = fWidth->Eval(MyGamma.gammaEn);
				for(int j = 0; j < counts; j++){
						hEscPeaks[source_number]->Fill( gRandom->Gaus(MyGamma.gammaEn-511.,width) );
				}
		}
}


//this function is used to produce the full simulated spectrum
//this function reads in the simulated Full energy peaks, single escape peaks and the background extracted from the real spectrum
//this function creates three new spectra
//One adds the full-energy peaks, single esc, peaks and the background, this is the full simulated spectrum
//Another adds the full energy peaks to the background
//The third adds the escape peaks to the background
void BuildSimuledSpectra(){


	const int Nbins = hBkgr->GetXaxis()->GetNbins();
	double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
	double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
	
	hFullSim = new TH1D("hFullSim","Full Sim Spectrum",Nbins, x_low, x_max);
	hFullSim->SetLineColor(kRed);
	hFullSim->Add(hBkgr);
	
	int hist_colors[] = {6, 417, 1, 900-4, 432, 801, 880, 861, 625, 416};
	for(int i = 0; i < used_sources; i++){
		cout << i << endl;
		hSimSource[i] = new TH1D( Form("hSim_%s", source_name.at(i).c_str() ),  Form("Simulated Source Decay on Background: %s",  source_name.at(i).c_str() ), Nbins, x_low, x_max);
		hSimSource[i]->SetLineColor( hist_colors[i] );
		hSimSource[i]->Add(hBkgr);
		hSimSource[i]->Add(hSimPeaks[i]);
		hSimSource[i]->Add(hEscPeaks[i]);
		hFullSim->Add(hSimPeaks[i]);
		hFullSim->Add(hEscPeaks[i]);
	}
	

	hRealSpectra->SetFillColor(kBlue);
	hRealSpectra->SetFillStyle(3003);		
	hRealSpectra->Draw("hist");
	hFullSim->Draw("histsame");
	for(int i = 0; i < used_sources; i++){
		hSimSource[i]->Draw("histsame");
	}	
	hBkgr->Draw("histsame");
}


void WriteSimulation(string filename){
		
		TFile *fNew = TFile::Open(filename.c_str(), "RECREATE");
		hRealSpectra->Write();
		hFullSim->Write();
		hBkgr->Write();
		for(int i = 0; i < used_sources; i++){
			hSimSource[i]->Write();
		}
}
