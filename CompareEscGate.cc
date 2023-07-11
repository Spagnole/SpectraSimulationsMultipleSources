#include "../ReadInputData.cc"

TH1D *hGatePeak; TH1D *hGateEspPeak;
TH1D *hGateBkgr; TH1D *hGateEscBkgr;
TH1D *hGateEspPeakScaled;

void Gate(TH2D *mat, int PeakLow, int PeakHigh, int BkgrLow, int BkgrHigh, int EscBkgrLow, int EscBkgrHigh){

	hGatePeak = mat->ProjectionX(Form("hPeak_%d_%d",PeakLow,PeakHigh),PeakLow,PeakHigh);
	hGateBkgr = mat->ProjectionX(Form("hBkgr_%d_%d",BkgrLow,BkgrHigh),BkgrLow,BkgrHigh);
	hGateEspPeak = mat->ProjectionX(Form("hEscPeak_%d_%d",PeakLow+511,PeakHigh+511),PeakLow+511,PeakHigh+511);
	hGateEscBkgr = mat->ProjectionX(Form("hEscBkgr_%d_%d",EscBkgrLow,EscBkgrHigh),EscBkgrLow,EscBkgrHigh);
	
	double ComptonPeak = hBkgr->Integral(PeakLow,PeakHigh);
	double ComptonBkgr = hBkgr->Integral(BkgrLow,BkgrHigh);
	double ComptonEsc = hBkgr->Integral(PeakLow+511,PeakHigh+511);
	double ComptonEscBkgr = hBkgr->Integral(EscBkgrLow,EscBkgrHigh);
	
	hGatePeak->SetTitle( Form("Gate: %d-%d",PeakLow,PeakHigh) );
	
	hGatePeak->Add(hGateBkgr,-1*ComptonPeak/ComptonBkgr );
	hGateEspPeak->Add(hGateEscBkgr,-1*ComptonEsc/ComptonEscBkgr );
	hGateEspPeak->SetLineColor(kRed);
	hGateEspPeakScaled = (TH1D*)hGateEspPeak->Clone();
	hGateEspPeakScaled->SetName( Form("hEscPeakScaled_%d_%d",PeakLow+511,PeakHigh+511) );
	hGateEspPeakScaled->SetTitle( Form("Gate (Scaled):%d_%d",PeakLow+511,PeakHigh+511) );
	hGateEspPeakScaled->SetLineColor(6);
	if( fEscPeak != NULL ){
		double scale_par = fEscPeak->Eval( 0.5*(PeakLow+511 + PeakHigh+511) );
		cout << "Scaling Histogram by: " << scale_par << endl;
		hGateEspPeakScaled->Scale( 1*scale_par );
	}
	else{
		cout << "No function for escape peak intensity given!\nNo scaling of single escape peak coincidence performed!\n";
		cout << "Did you execute ReadEscapePeaks() function?\n";
				
	}
	hGatePeak->Draw("hist");
	hGateEspPeakScaled->Draw("histsame");
//	new TCanvas(); hGatePeak->Draw("hist");
//	new TCanvas(); hGateBkgr->Draw("hist");
//	new TCanvas(); hGateEspPeak->Draw("hist");
//	new TCanvas(); hGateEscBkgr->Draw("hist");

}



