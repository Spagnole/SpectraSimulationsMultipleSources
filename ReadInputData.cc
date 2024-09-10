#include "ReadInputData.hh"


//////////////////////////////////////////////////////////////
///// PLEASE NOT ALL ENERGIES SHOULD BE GIVEN IN keV!!! //////
//////////////////////////////////////////////////////////////

//This is for the list of transitions from your decay scheme
//there are six components to this vector
//Initial level energy ---> gamma-ray energy ---> final level energy ---> g.-ray intensity ---> g-ray int. uncertainty ---> level population int. ---> g.-ray Branching ratio
//see comments for ReadDecayScheme() for more details


//this function calls for your input decay scheme
//This function requires a 5 column text file as the input
//The input is as follows
//Initial level energy ---> gamma-ray energy ---> final level energy ---> gamma-ray intensity ---> gamma-ray intensity uncertainty
//the gamma rays are added to the vector 'AssignedTransition'
//the 'level population int.' and 'g.-ray Branching ratio' components are only used for simulations of coincidence spectra
//these values are not read in by this function but are calculated by the function CalcLevelFeedingAndGammaBR()
//this function is included in the SimulateCoincidences.cc code

void ReadDecayScheme(string filename = "example-152Eu/152Gd.dat", string enter_source_name = Form("source_%d",used_sources)){
		ifstream infile( filename.c_str() );
		if( !infile ){
				cout << filename << " is not open!" << endl;
				return;
		}
		int counter = 0 ;
		string line;
		while (getline(infile, line)) {
				// Skip lines starting with '#'
				if (line.empty() || line[0] == '#') {
						cout << "Skipping line " << line << endl;
						continue;
				}
				
				stringstream ss(line);
				double a[5] = {0}; string comment = "";
				MyTransition MyGamma;
				ss >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> comment;
				MyGamma.Index = counter;
				MyGamma.lvlEn = a[0];
				MyGamma.gammaEn = a[1];
				MyGamma.finalLvl = a[2];
				MyGamma.gammaInt = a[3];
				MyGamma.gammaIntError = a[4];
				MyGamma.lvlPop = 0;
				MyGamma.gammaBr = 0;
				gammaList[used_sources].push_back(MyGamma);
				counter++;
		}
		used_sources++;
		source_name.push_back(enter_source_name);
		infile.close();
}



		


//The GetEfficiency() function is used to read in your gamma-ray efficiency and add produce a graph of efficiency as a function of energy
//this function reads in a two column text file
//the first column is energy in keV
//the second column is the efficiency at that energy
//This efficiency graph is used to correct intensities when filling the simulated spectra
void GetEfficiency(string eff_filename = "MyExpEffnew.dat"){
		ifstream myfitresult( eff_filename.c_str() );
		if( !myfitresult.is_open() ){
				cout << eff_filename << " is not open!" << endl;
				return;
		}
		
		gEff = new TGraph(); //this is the graph for the efficiency curve
		
		string line;
		while (getline(myfitresult, line)) {
				// Skip lines starting with '#'
				if (line.empty() || line[0] == '#') {
						cout << "Skipping line " << line << endl;
						continue;
				}
				stringstream ss(line);
				double a[2] = {0}; string comments = "";
				ss >> a[0] >> a[1] >> comments;
				gEff->SetPoint( gEff->GetN(), a[0], a[1]);
		}
		myfitresult.close();
}



//the GetPeakWidth() function is used to set the peak width as a function of energy for your simulated spectra
//This file requires a four column text file, the input should be the following
//gamma-ray energy ---> g.-ray energy error ---> peak width ---> peak width error
//please use term sigma! NOT FWHM!!!!
//this function fills a graph with peak widths as a function of energy
//a linear function is then fit to the data
//this function is used to get the peak widths for the simulated spectra
void GetPeakWidth(string peak_widths_filename = "PeakWidths.dat"){

		ifstream input( peak_widths_filename.c_str() );
		if( !input.is_open() ){
			cout << peak_widths_filename << " is not open!" << endl;
			return;
		}
		
		gSigma = new TGraphErrors();	gSigma->SetName("gSigma");	gSigma->SetMarkerStyle(20);	gSigma->SetMarkerColor(kBlue);

		string line;
		while (getline(input, line)) {
				// Skip lines starting with '#'
				if (line.empty() || line[0] == '#') {
						cout << "Skipping line " << line << endl;
						continue;
				}
				stringstream ss(line);
				double a[4] = {0}; string comments = "";
				ss >> a[0] >> a[1] >> a[2] >> a[3] >> comments;
				gSigma->SetPoint( gSigma->GetN(), a[0], a[2]);
				gSigma->SetPointError( gSigma->GetN()-1,  a[1], a[3]);
		}

		fWidth = new TF1("fWidth","[0]+[1]*x",0,8000);
		fWidth->SetParameters(9.57477e-01, 2.59267e-04);
		gSigma->Fit("fWidth");
		//gSigma->Draw("AP");
		//fWidth->Draw("same");
		input.close();
}

void DefaultPeakWidths(double offset = 9.57477e-01, double gain = 2.59267e-04){
		fWidth = new TF1("fWidth","[0]+[1]*x",0,8000);
		fWidth->SetParameters(offset,gain);
}

//the ReadEscapePeaks() function is used to get escape peak intensites relative to the full energy peak as a function of energy
//This file requires a four column text file, the input should be the following
//gamma-ray energy ---> g.-ray energy error ---> Esc.-Peak Int. / Full-Energy Peak Int  ---> uncertainty (Esc.-Peak Int. / Full-Energy Peak Int)
//this function is required for the simulation of escape peak intensities
void ReadEscapePeaks(string EscPeaksFilename = "EscapePeaks.dat"){
		
		gEscPeaks = new TGraphErrors();	gEscPeaks->SetName("gSigma");	gEscPeaks->SetMarkerStyle(20);	gEscPeaks->SetMarkerColor(kBlue);
		
		ifstream input( EscPeaksFilename.c_str() );
		if( !input.is_open() ){
				cout << EscPeaksFilename << " is not open!" << endl;
				return;
		}
		
		string line;
		while (getline(input, line)) {
				// Skip lines starting with '#'
				if (line.empty() || line[0] == '#') {
						cout << "Skipping line " << line << endl;
						continue;
				}
				stringstream ss(line);
				double a[4] = {0}; string comments = "";
				ss >> a[0] >> a[1] >> a[2] >> a[3] >> comments;
				gEscPeaks->SetPoint( gEscPeaks->GetN(), a[0], a[2]);
				gEscPeaks->SetPointError( gEscPeaks->GetN()-1,  a[1], a[3]);
				if(PrintEscPeakData) cout << line << endl;
		}
			
		fEscPeak = new TF1("fEscPeak","[0]+[1]*x+[2]*x*x",0,8000);
		fEscPeak->SetParameters(-1.45859e-02, 1.11648e-06, 7.51546e-09);
		gEscPeaks->Fit("fEscPeak");
		//gEscPeaks->Draw("AP");
		//fEscPeak->Draw("same");
		input.close();
}

void DefaultEscapePeaks(double offset = -1.45859e-02, double gain = 1.11648e-06, double quadratic = 7.51546e-09){
		fEscPeak = new TF1("fEscPeak","[0]+[1]*x+[2]*x*x",0,8000);
		fEscPeak->SetParameters(offset, gain, quadratic);
}

//this function is used to extract a background from your experimntal spectrum
//the default values of the function will extract a rather crude background
//the user is free to play with these parameters to try and improve the extracted background
void GetSpectrumBackground(TH1D *h, int iterations = 50, int decreasewindow = 1, int backorder = 2, bool smoothing = false, int smoothwindow = 3, bool compton = 0){

	//extract binning information from the experimental spectrum
	//we want to use the same binning to compare
	const int Nbins = h->GetXaxis()->GetNbins();
	double x_low = h->GetXaxis()->GetBinLowEdge(1);
	double x_max = h->GetXaxis()->GetBinUpEdge(Nbins);
	double source[Nbins];
	
	hBkgr = new TH1D("hBkgr","Simulated Background",Nbins,x_low,x_max);
	
	TSpectrum *s = new TSpectrum();
	for (int i = 0; i < Nbins; i++) source[i]=h->GetBinContent(i + 1);
	//s->Background(source,Nbins,75,TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2,kFALSE, TSpectrum::kBackSmoothing3,kFALSE);
	s->Background(source,Nbins,iterations,decreasewindow,backorder,smoothing,smoothwindow,compton);
	for (int i = 0; i < Nbins; i++) hBkgr->SetBinContent(i + 1,source[i]);      
	h->Draw("hist");
	hBkgr->SetLineColor(kOrange+1);
	hBkgr->Draw("SAME L");
                 
}

//this function is used to get your real experimental spectra
//user should provide the name of the root file and the name of the spectrum
void GetRealSpectra(string rootfilename = "ExampleFile.root", string RealHistName = "hgE_56Co"){

	//open root file
	if(gSystem->AccessPathName(rootfilename.c_str())){ //check to see root file exists!
		cout << "The file " << rootfilename << " doesn't exists\n";
		return;
	}
	else fRealData = TFile::Open( rootfilename.c_str() ); //opening root file
	
	if( !fRealData->GetListOfKeys()->Contains( RealHistName.c_str() )){ //check to see histogram exists!
		cout << "The histogram " << RealHistName << " doesn't exists\n";
		return;
	}	
	else{ //getting root spectrum
		hRealSpectra = (TH1D*)fRealData->Get( RealHistName.c_str() );
		hRealSpectra->Draw("hist");
		hRealSpectra->GetXaxis()->SetTitle("Energy (keV)");
		hRealSpectra->GetXaxis()->CenterTitle();
	}
}

//Addition info for GetSpectrumBackground() function
/*
       kBackOrder2 =0,
       kBackOrder4 =1,
       kBackOrder6 =2,
       kBackOrder8 =3,
       kBackIncreasingWindow =0,
       kBackDecreasingWindow =1,
       kBackSmoothing3 =3,
       kBackSmoothing5 =5,
       kBackSmoothing7 =7,
       kBackSmoothing9 =9,
       kBackSmoothing11 =11,
       kBackSmoothing13 =13,
       kBackSmoothing15 =15
       
       Background(Double_t* spectrum, Int_t ssize, Int_t numberIterations, Int_t direction, Int_t filterOrder, bool smoothing, Int_t smoothWindow, bool compton)
*/

void printGammaList(int SourceNum = 0){
		for(auto MyGamma : gammaList[SourceNum] ) MyGamma.Display();
}

void printLevels(int source_number = 0){
		for(auto level : levelList[source_number]) level.Display();
}
void printLevelsReverse(int source_number = 0){
		for(auto level = levelList[source_number].rbegin(); level != levelList[source_number].rend(); level++)
				level->Display();
}

void printCoincidences(int source_number){
		for(auto gg : coincList[source_number]) gg.Display();
}

/*
void OldReadDecayScheme(string filename = "example-152Eu/152Gd.dat", string enter_source_name = Form("source_%d",used_sources)){
		ifstream input( filename.c_str() );
		if( !input.is_open() ){
				cout << filename << " is not open!" << endl;
				return;
		}
		double a[5];
		MyTransition MyGamma;
		input >> a[0] >> a[1] >> a[2] >> a[3] >> a[4];
		while( !input.eof() ){
				MyGamma.lvlEn = a[0];
				MyGamma.gammaEn = a[1];
				MyGamma.finalLvl = a[2];
				MyGamma.gammaInt = a[3];
				MyGamma.gammaIntError = a[4];
				MyGamma.lvlPop = 0;
				MyGamma.gammaBr = 0;
				gammaList[used_sources].push_back(MyGamma);
				if(PrintReadData) cout << a[0] << "\t" <<  a[1] << "\t" <<
						a[2] << "\t" <<  a[3] << "\t" <<  a[4] << endl;
				input >> a[0] >> a[1] >> a[2] >> a[3] >> a[4];
		}
		used_sources++;
		source_name.push_back(enter_source_name);
}
*/
/*
void OldGetEfficiency(string eff_filename = "MyExpEffnew.dat"){
		ifstream myfitresult( eff_filename.c_str() );
		if( !myfitresult.is_open() ){
				cout << eff_filename << " is not open!" << endl;
				return;
		}
		gEff = new TGraph();
		double a[2];
		myfitresult >> a[0] >> a[1];
		while( !myfitresult.eof() ){
				gEff->SetPoint( gEff->GetN(), a[0], a[1]);
				myfitresult >> a[0] >> a[1];
		}
}
*/
/*
void OldGetPeakWidth(string peak_widths_filename = "PeakWidths.dat"){

	gSigma = new TGraphErrors();
	gSigma->SetName("gSigma");
	gSigma->SetMarkerStyle(20);
	gSigma->SetMarkerColor(kBlue);
	
	ifstream input( peak_widths_filename.c_str() );
	if( !input.is_open() ){
		cout << peak_widths_filename << " is not open!" << endl;
		return;
	}
	double a[4];
	input >> a[0] >> a[1] >> a[2] >> a[3];
	while( !input.eof() ){
		if( a[2] > a[3] ){
			gSigma->SetPoint( gSigma->GetN(), a[0], a[2]);
			gSigma->SetPointError( gSigma->GetN()-1,  a[1], a[3]);
			if(PrintPeakWidths) cout << a[0] << "\t" << a[1] << "\t" << a[2] << "\t" << a[3] << endl;
		}
			input >> a[0] >> a[1] >> a[2] >> a[3];
		
	}
	fWidth = new TF1("fWidth","[0]+[1]*x",0,8000);
	fWidth->SetParameters(9.57477e-01, 2.59267e-04);
	//gSigma->Draw("AP");
	gSigma->Fit("fWidth");
	//fWidth->Draw("same");
}
*/
