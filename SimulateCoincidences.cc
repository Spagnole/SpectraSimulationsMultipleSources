#include "ReadInputData.cc"

bool PrintReducedCoincList = false;

//this vector is a list of levels, this vector is populated by the GetLevelList() function.
//The components of this vector are: Level Index ---> Level Energy ---> Level Population
vector<tuple<int, double, double>> LevelEnergy[NSources];

//list of gamma-gamma coincidence intensities
//The components of this vector are: Gamma 1 ---> Gamma 2 ---> Intensity of coincidence
vector<tuple<double,double,double>> gg_coinc[NSources]; 

TH2D *sim_gg_mat[NSources];
TH1D *SimMatProj[NSources];
TH1D *SimEscMatProj[NSources];
TH1D *SimEscGatedSpectra[NSources];
//TH1D *hEscPeaks[NSources];
//TH1D *hSourceSimSpectrum[NSources];
//TH1D *hFullSim;


//this function gets a list of levels from the 'AssignedTransition' vector from the ReadInputData.cc file
//this function requires that ReadDecayScheme() function has been utilized
void GetLevelList(int source_number){

	//for(int source_number = 0; source_number < used_sources; source_number++){

		if( AssignedTransition[source_number].size() == 0 ){
			cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
			return;
		}

		double prev_energy = -1;
		int counter = 0;
		for(int i = 0; i < AssignedTransition[source_number].size(); i++){
			if( get<0>(AssignedTransition[source_number].at(i)) != prev_energy ){
				LevelEnergy[source_number].push_back( make_tuple(counter, get<0>(AssignedTransition[source_number].at(i)), 0) );
				counter++;
			}
			prev_energy = get<0>(AssignedTransition[source_number].at(i));
		}
		
		for(int i = 0; i < LevelEnergy[source_number].size(); i++){
			if(PrintLevelList) cout << get<0>(LevelEnergy[source_number].at(i)) << "\t" << get<1>(LevelEnergy[source_number].at(i)) <<endl;
		}
	//}
}

//this function is used to correct the 'AssignedTransition' list
//this correction is required if your final level energies do not exactly match the level energies in the 'LevelEnergy' vector

void FixTransitionList(int source_number){

	//for(int source_number = 0; source_number < used_sources; source_number++){
		if( AssignedTransition[source_number].size() == 0 ){
			cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
			return;
		}

		double min_diff;
		double difference;
		double new_level_energy;
		for(int i = 0; i <  AssignedTransition[source_number].size(); i++){
			min_diff = 10000;
			for(int j = 0; j < LevelEnergy[source_number].size(); j++){
				difference = TMath::Abs( get<2>(AssignedTransition[source_number].at(i)) - get<1>(LevelEnergy[source_number].at(j)) );
				if( difference < min_diff ){
					new_level_energy = get<1>(LevelEnergy[source_number].at(j));
					min_diff = difference;
				}
			}
			get<2>(AssignedTransition[source_number].at(i)) = new_level_energy;		
		}
	//}
}


//this function calculates the level population of each state and the gamma-ray branching ratio of each transition
//this is function is required to properly calculate the expected gamma-gamma coincidence intensity
void CalcLevelFeedingAndGammaBR(int source_number){
	
	//for(int source_number = 0; source_number < used_sources; source_number++){
		double sum_temp;
		for(int i = 0; i < LevelEnergy[source_number].size(); i++){
			sum_temp = 0;
			for(int j = 0; j <  AssignedTransition[source_number].size(); j++){
				if(  get<1>(LevelEnergy[source_number].at(i)) == get<0>(AssignedTransition[source_number].at(j))){
					sum_temp += get<3>(AssignedTransition[source_number].at(j));
				}
			}
			get<2>(LevelEnergy[source_number].at(i)) = sum_temp;
		}
		
		//calculating gamma-ray branching ratios
		for(int i = 0; i <  AssignedTransition[source_number].size(); i++){
			if( i == 0 ) continue;
			for(int j = 0; j < LevelEnergy[source_number].size(); j++){
				if(  get<1>(LevelEnergy[source_number].at(j)) == get<0>(AssignedTransition[source_number].at(i))){
					get<4>(AssignedTransition[source_number].at(i)) = get<2>(LevelEnergy[source_number].at(j));
					get<5>(AssignedTransition[source_number].at(i)) = get<3>(AssignedTransition[source_number].at(i)) / get<2>(LevelEnergy[source_number].at(j));
				}
			}
		}
		
		//calculate level population
		for(int i = 0; i <  AssignedTransition[source_number].size(); i++){
			if( i == 0 ) continue;		
			for(int j = 0; j < AssignedTransition[source_number].size(); j++){
				if(  get<2>(AssignedTransition[source_number].at(j)) == get<0>(AssignedTransition[source_number].at(i))){
					get<4>(AssignedTransition[source_number].at(i)) = get<4>(AssignedTransition[source_number].at(i)) - get<3>(AssignedTransition[source_number].at(j));
				}
			}
		}
		
		//calculating gamma-ray branching ratios
		for(int i = AssignedTransition[source_number].size()-1; i >= 0; i--){
			if(PrintCalcBR) cout << get<0>(AssignedTransition[source_number].at(i)) << "\t" << get<1>(AssignedTransition[source_number].at(i)) << "\t" << get<2>(AssignedTransition[source_number].at(i)) << "\t";
			if(PrintCalcBR) cout << get<3>(AssignedTransition[source_number].at(i)) << "\t" << get<4>(AssignedTransition[source_number].at(i)) << "\t" << get<5>(AssignedTransition[source_number].at(i)) << endl;
		}
	//}
}


void AddToCoincidenceList(int source_number, vector<int> feeder_index, int decay_index, double intensity ){

	double decay_energy = get<1>(AssignedTransition[source_number].at(decay_index));
	for(int i = 0; i < feeder_index.size(); i++){
		double feeder_energy = get<1>(AssignedTransition[source_number].at( feeder_index.at(i) ) );	
		gg_coinc[source_number].push_back( make_tuple(feeder_energy,decay_energy, intensity) );
	}
}

void PrintUnfinishedCascade(int source_number, vector<int> list){

	cout  << "WARNING!!! Cascade not finished"<<endl;	
	for(int x = 0; x < list.size(); x++){
		int y = list.at(x);
		for(int z = 0; z < x; z++) cout <<"\t";
		cout << get<0>(AssignedTransition[source_number].at(y)) << " --> " << 
		get<1>(AssignedTransition[source_number].at(y)) << " --> " << get<2>(AssignedTransition[source_number].at(y)) << endl;
	}
	cout << "User can add nested loops for FindCoincidences() function!\n";
}

//this function calculates the expected number of gamma-gamma coincidences
//this function uses a large number of nested loops (Should be replaced by some recursive loop?)
//if a gamma-gamma coincidence is separated by about 5 intermediate gamma-rays it will not be added to the list of coincidences
//if you require such coincidences ---> Add more nested loops
void FindCoincidences(int source_number){
		
	if( AssignedTransition[source_number].size() == 0 ){
		cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
		return;
	}

	double NDecays, Ngg[10];
	
	gg_coinc[source_number].clear();
	
	for(int i = AssignedTransition[source_number].size()-1; i >= 0 ; i--){
		if( get<2>(AssignedTransition[source_number].at(i))  == 0 ) continue;
		NDecays = get<4>(AssignedTransition[source_number].at(i))*get<5>(AssignedTransition[source_number].at(i));
		for(int j = i-1; j >= 0 ; j--){
			if( get<2>(AssignedTransition[source_number].at(i)) == get<0>(AssignedTransition[source_number].at(j))){
				Ngg[0] = NDecays*get<5>(AssignedTransition[source_number].at(j));
				AddToCoincidenceList(source_number,{i},j, Ngg[0] );
				if( get<2>(AssignedTransition[source_number].at(j))  == 0 ) continue;
				for(int k = j-1; k >= 0 ; k--){
					if( get<2>(AssignedTransition[source_number].at(j)) == get<0>(AssignedTransition[source_number].at(k))){
						Ngg[1] = Ngg[0]*get<5>(AssignedTransition[source_number].at(k));
						AddToCoincidenceList(source_number,{i,j},k, Ngg[1] );
						if( get<2>(AssignedTransition[source_number].at(k))  == 0 ) continue;
						for(int l = k-1; l >= 0 ; l--){
							if( get<2>(AssignedTransition[source_number].at(k)) == get<0>(AssignedTransition[source_number].at(l))){
								Ngg[2] = Ngg[1]*get<5>(AssignedTransition[source_number].at(l));
								AddToCoincidenceList(source_number,{i,j,k},l, Ngg[2] );
								if( get<2>(AssignedTransition[source_number].at(l))  == 0 ) continue;
								for(int m = l-1; m >= 0 ; m--){
									if( get<2>(AssignedTransition[source_number].at(l)) == get<0>(AssignedTransition[source_number].at(m))){
										Ngg[3] = Ngg[2]*get<5>(AssignedTransition[source_number].at(m));
										AddToCoincidenceList(source_number,{i,j,k,l},m, Ngg[3] );
										if( get<2>(AssignedTransition[source_number].at(m))  == 0 ) continue;
										for(int n = m-1; n >= 0 ; n--){
											if( get<2>(AssignedTransition[source_number].at(m)) == get<0>(AssignedTransition[source_number].at(n))){
												Ngg[4] = Ngg[3]*get<5>(AssignedTransition[source_number].at(n));
												AddToCoincidenceList(source_number,{i,j,k,l,m},n, Ngg[4] );
												if( get<2>(AssignedTransition[source_number].at(n))  == 0 ) continue;
												for(int p = n-1; p >= 0 ; p--){
													if( get<2>(AssignedTransition[source_number].at(n)) == get<0>(AssignedTransition[source_number].at(p))){
														Ngg[5] = Ngg[4]*get<5>(AssignedTransition[source_number].at(p));
														AddToCoincidenceList(source_number,{i,j,k,l,m,n},p, Ngg[5] );
														if( get<2>(AssignedTransition[source_number].at(p)) !=0) PrintUnfinishedCascade(source_number,{i,j,k,l,m,n});
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


//this function is used to reduce the size of the gamma-gamma coincidence list
//if a given coincidence is given a number of times, this function will reduce the entries of that coincidence to one instance but adds all of the intensities
//strictly speaking this function is not absolutely required for the tool to function! BUT YOU SHOULD STILL USE IT!
void ReduceGammaGammaList(int source_number){
	
	for(int i = 0; i < gg_coinc[source_number].size();i++){
		for(int j = 0; j < gg_coinc[source_number].size(); j++){
			if(i==j) continue;
			if( get<0>(gg_coinc[source_number].at(i)) == get<0>(gg_coinc[source_number].at(j)) && get<1>(gg_coinc[source_number].at(i)) == get<1>(gg_coinc[source_number].at(j)) ){
				if( PrintReducedCoincList ) cout << get<0>(gg_coinc[source_number].at(i)) << "\t" << get<1>(gg_coinc[source_number].at(i)) << "\t"  << get<2>(gg_coinc[source_number].at(i)) <<endl;
				if( PrintReducedCoincList ) cout <<"\t" << get<0>(gg_coinc[source_number].at(j)) << "\t" << get<1>(gg_coinc[source_number].at(j)) << "\t"  << get<2>(gg_coinc[source_number].at(j)) <<endl;
				get<2>(gg_coinc[source_number].at(i)) =  get<2>(gg_coinc[source_number].at(i)) + get<2>(gg_coinc[source_number].at(j));
				if( PrintReducedCoincList ) cout << "\t\t" << get<0>(gg_coinc[source_number].at(i)) << "\t" << get<1>(gg_coinc[source_number].at(i)) << "\t"  << get<2>(gg_coinc[source_number].at(i)) <<endl;
				gg_coinc[source_number].erase(gg_coinc[source_number].begin()+j);
				j=j-1;
			}
		}
	}
}

void PrintCoincidences(int source_number, string filename){
	
	ofstream outfile( filename.c_str() );
	for(int i = 0; i < gg_coinc[source_number].size();i++){
		outfile << get<0>(gg_coinc[source_number].at(i)) << "\t" << get<1>(gg_coinc[source_number].at(i)) << "\t"  << get<2>(gg_coinc[source_number].at(i)) <<endl;
	}
}

//this function fills a simulated gamma-gamma coincidence matrix
//user can provide specified binning of the matrix
//this function requires that the list of coincidences has been produced
//as well as reading in the gamma-ray efficiency and peak widths
//
void FillCoincMatrix(int source_number, int NBins = -1, double low = -1, double upp = -1){
	
	if( NBins == -1){
		NBins = hRealSpectra->GetXaxis()->GetNbins();
		low = hRealSpectra->GetXaxis()->GetBinLowEdge(1);
		upp = hRealSpectra->GetXaxis()->GetBinUpEdge(NBins);
	}
	//Form("hPeaks_%s", source_name.at(source_number).c_str()),Form("Simulated Source Peaks: %s",  source_name.at(source_number).c_str() )
	
	sim_gg_mat[source_number] = new TH2D(Form("sim_gg_mat_%s", source_name.at(source_number).c_str() ), Form("Simulated g-g matrix Source: %s",source_name.at(source_number).c_str() ),NBins,low,upp,NBins,low,upp);
	double gamma1,gamma2,sigma1,sigma2;
	double counts;
	for(int i = 0; i < gg_coinc[source_number].size(); i++){
		gamma1 = get<0>(gg_coinc[source_number].at(i));
		gamma2 = get<1>(gg_coinc[source_number].at(i));
		counts = (14./15) * get<2>(gg_coinc[source_number].at(i)) * gEff->Eval(gamma1) * gEff->Eval(gamma2); //user may require different normalisation if using relative efficiency curve
		sigma1 = fWidth->Eval(gamma1);
		sigma2 = fWidth->Eval(gamma2);
		for(int j = 0; j < counts; j++){
			sim_gg_mat[source_number]->Fill( gRandom->Gaus(gamma1,sigma1), gRandom->Gaus(gamma2,sigma2));
			sim_gg_mat[source_number]->Fill( gRandom->Gaus(gamma2,sigma2), gRandom->Gaus(gamma1,sigma1));
		}
	}
	
}


//this function gates on the simulated gamma-gamma matrix to produce a coincidence spectrum
//used should use the same gate they used for their experimental efficiency curve
void GateOnSimMat(int source_number, int low, int up){
	
	hSimPeaks[source_number]  = sim_gg_mat[source_number]->ProjectionY(Form("SimProj_%s_%d_%d",source_name.at(source_number).c_str(),low,up),low,up);
	hSimPeaks[source_number]->Draw("hist");
	hSimPeaks[source_number]->SetLineColor(kRed);
	hSimPeaks[source_number]->SetFillColor(kRed);
	hSimPeaks[source_number]->SetFillStyle(3000);
	
	double MeanEnergy = ( (double)low + (double)up )/2.;
	if( MeanEnergy > 1000. && gEscPeaks != NULL){
		SimEscMatProj[source_number]  = sim_gg_mat[source_number]->ProjectionY(Form("SimEscProj_%s_%d_%d",source_name.at(source_number).c_str(),low+511,up+511),low+511,up+511);		
		SimEscMatProj[source_number]->Scale(fEscPeak->Eval(MeanEnergy));
		//SimEscGatedSpectra[source_number] = SimEscMatProj[source_number]->Clone();
		SimEscGatedSpectra[source_number] = new TH1D(*SimEscMatProj[source_number]);
		SimEscGatedSpectra[source_number]->SetName( Form("SimGateOnEscPeak_%s_%d_%d",source_name.at(source_number).c_str(),low+511,up+511) );
		SimEscGatedSpectra[source_number]->SetTitle( Form("Sim. Gated Spectra Gate on Escape Peak Source:%s, %d - %d",source_name.at(source_number).c_str(),low+511,up+511) );
	}
}

//this function fills the simulated escape peaks 
//this is done by using the output TH1D from the GateOnSimMat() function
void FillEscPeakSpec(int source_number){

	if( gEscPeaks == NULL ) return;

	int Nbins = hSimPeaks[source_number]->GetXaxis()->GetNbins();
	double low = hSimPeaks[source_number]->GetXaxis()->GetBinLowEdge(1);
	double upp = hSimPeaks[source_number]->GetXaxis()->GetBinUpEdge(Nbins);
	
	hEscPeaks[source_number] = new TH1D(Form("hEscPeaks_%s", source_name.at(source_number).c_str()),Form("Simulated Source Single-Escape Peaks: %s",  source_name.at(source_number).c_str() ),Nbins, low, upp);
	hEscPeaks[source_number]->SetLineColor(6);
	hEscPeaks[source_number]->SetFillColor(6);
	hEscPeaks[source_number]->SetFillStyle(3000);
	double gamma1,sigma1;
	double counts, scale; 
	
	for(int i = 1; i <= Nbins; i++){
		gamma1 = hSimPeaks[source_number]->GetBinCenter(i);
		if(gamma1 < 1500. ) continue;
		counts = hSimPeaks[source_number]->GetBinContent(i);
		for(int j = 0; j < counts; j++){
			hEscPeaks[source_number]->Fill( gamma1-511. , fEscPeak->Eval(gamma1) );
		}
	}
}

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
		if( hEscPeaks[i] != NULL ) hSimSource[i]->Add(hEscPeaks[i]);
		if( SimEscGatedSpectra[i] != NULL){
			SimEscGatedSpectra[i]->Add(hBkgr);
			SimEscGatedSpectra[i]->SetLineColor( hist_colors[i] );
			SimEscGatedSpectra[i]->SetFillColor( hist_colors[i] );
			SimEscGatedSpectra[i]->SetFillStyle(3144);
			hFullSim->Add( SimEscMatProj[i] );
		}		
		hFullSim->Add(hSimPeaks[i]);
		hFullSim->Add(hEscPeaks[i]);
	}
	

	hRealSpectra->SetFillColor(kBlue);
	hRealSpectra->SetFillStyle(3003);		
	hRealSpectra->Draw("hist");
	
	for(int i = 0; i < used_sources; i++){
		hSimSource[i]->Draw("histsame");
		SimEscGatedSpectra[i]->Draw("histsame");
	}
	hFullSim->Draw("histsame");	
	hBkgr->Draw("histsame");
}
/*
//this function is identical to the function used by SimulateSingles.cc (Perhaps it should be put in the ReadInputData.cc file?)
//this function is used to produce the full simulated spectrum
//this function reads in the simulated Full energy peaks, single escape peaks and the background extracted from the real spectrum
//this function creates three new spectra
//One adds the full-energy peaks, single esc, peaks and the background, this is the full simulated spectrum
//Another adds the full energy peaks to the background
//The third adds the escape peaks to the background
void BuildSimuledSpectra(TH1D *hSimPeaks, TH1D *hBkgr, TH1D *hEscPeaks){

	const int Nbins = hBkgr->GetXaxis()->GetNbins();
	double x_low = hBkgr->GetXaxis()->GetBinLowEdge(1);
	double x_max = hBkgr->GetXaxis()->GetBinUpEdge(Nbins);
	hFullSim = new TH1D("hFullSim","Full Sim Spectrum",Nbins, x_low, x_max);
	hSimPeaks_Bkgr = new TH1D("hSimPeaks_Bkgr","Sim Peaks on Bkgr",Nbins, x_low, x_max);
	hEscPeaks_Bkgr = new TH1D("hEscPeaks_Bkgr","Sim Esc, Peaks on Bkgr",Nbins, x_low, x_max);


	hFullSim->SetLineColor(kRed);
	hFullSim->Add(hBkgr);
	hFullSim->Add(hSimPeaks);
	hFullSim->Add(hEscPeaks);

	hSimPeaks_Bkgr->SetLineColor(kGreen+2);
	hSimPeaks_Bkgr->Add(hBkgr);
	hSimPeaks_Bkgr->Add(hSimPeaks);
		
	hEscPeaks_Bkgr->SetLineColor(6);
	hEscPeaks_Bkgr->Add(hBkgr);
	hEscPeaks_Bkgr->Add(hEscPeaks);

	hRealSpectra->SetFillColor(kBlue);
	hRealSpectra->SetFillStyle(3003);
		
	hRealSpectra->Draw("hist");
	hSimPeaks_Bkgr->Draw("histsame");
	hEscPeaks_Bkgr->Draw("histsame");		
	hFullSim->Draw("histsame");
	hBkgr->Draw("histsame");
}

*/



