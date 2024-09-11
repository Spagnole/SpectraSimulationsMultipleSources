#include "ReadInputData.cc"
#include "WriteData.cc"
#include "BuildLevelScheme.cc"

TH2D *sim_gg_mat[NSources];
TH1D *SimMatProj[NSources];
TH1D *SimEscMatProj[NSources];
TH1D *SimEscGatedSpectra[NSources];

//ofstream outfile("nested-output.dat");
void AddToCoincidenceList(int source_number, vector<MyTransition> feederList, MyTransition Decay, double intensity ){
		
		gammagamma gg;
		for( auto Feeder : feederList){
				gg =  {Decay.gammaEn, Feeder.gammaEn, intensity};
				//outfile << Decay.gammaEn << "\t" << Feeder.gammaEn << "\t" << intensity <<endl;
				coincList[source_number].push_back(gg);
		}
}

void FindCoincidenceRec(int source_number, vector<MyTransition> &feederList, MyTransition currentTransition, double currentIntensity, int index) {
		// Loop through previous transitions starting from the current index - 1
		if( currentTransition.finalLvl == 0 ) return;
		
		for (int i = index - 1; i >= 0; --i) {
				MyTransition nextTransition = gammaList[source_number].at(i);
				// If the final level of the current transition matches the energy level of the next transition
				if (currentTransition.finalLvl == nextTransition.lvlEn) {
						double nextIntensity = currentIntensity * nextTransition.gammaBr;
						// Recursive call: Add the current transition to the feeder list and continue
						feederList.push_back(currentTransition);
						// Add to coincidence list
						AddToCoincidenceList(source_number, feederList, nextTransition, nextIntensity);
						FindCoincidenceRec(source_number, feederList, nextTransition, nextIntensity, i);
						feederList.pop_back(); // Backtrack
				}
		}
}


void FindCoincidences(int source_number) {
		if (gammaList[source_number].size() == 0) {
				cout << "No data in decay list!\nHave you read in the list of gamma-rays?" << endl;
				return;
		}
		coincList[source_number].clear();  // Clear the coincidence list for this source
		// Iterate through each transition as the starting point
		for (int g0 = gammaList[source_number].size() - 1; g0 >= 0; --g0) {
				MyTransition firstTransition = gammaList[source_number].at(g0);
				if (firstTransition.finalLvl == 0) continue;  // Skip transitions ending at level 0
				double NDecays = firstTransition.lvlPop * firstTransition.gammaBr;
				// Initialize an empty feeder list
				vector<MyTransition> feederList;
				// Start recursive search with current transition as the initial feeder
				FindCoincidenceRec(source_number, feederList, firstTransition, NDecays, g0);
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
		for(auto gg : coincList[source_number]){
				counts = ( (NDet-1) /NDet) * gg.coincidences * gEff->Eval(gg.gammaOne) * gEff->Eval(gg.gammaTwo);
				sigma1 = fWidth->Eval(gg.gammaOne);
				sigma2 = fWidth->Eval(gg.gammaTwo);
				for(int j = 0; j < counts; j++){
					sim_gg_mat[source_number]->Fill( gRandom->Gaus(gg.gammaOne,sigma1), gRandom->Gaus(gg.gammaTwo,sigma2));
					sim_gg_mat[source_number]->Fill( gRandom->Gaus(gg.gammaTwo,sigma2), gRandom->Gaus(gg.gammaOne,sigma1));
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
			if( SimEscGatedSpectra[i] != NULL) SimEscGatedSpectra[i]->Draw("histsame");
	}
	hFullSim->Draw("histsame");	
	hBkgr->Draw("histsame");
}

void WriteSimulation(string filename){
		
		TFile *fNew = TFile::Open(filename.c_str(), "RECREATE");
		hRealSpectra->Write();
		hFullSim->Write();
		hBkgr->Write();
		for(int i = 0; i < used_sources; i++){
			hSimSource[i]->Write();
				if( SimEscGatedSpectra[i] != NULL)
						SimEscGatedSpectra[i]->Write();
		}
}


//this function calculates the expected number of gamma-gamma coincidences
//this function uses a large number of nested loops (Should be replaced by some recursive loop?)
//if a gamma-gamma coincidence is separated by about 5 intermediate gamma-rays it will not be added to the list of coincidences
//if you require such coincidences ---> Add more nested loops
/*void FindCoincidences(int source_number){
		
		if( gammaList[source_number].size() == 0 ){
			cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
			return;
		}

		double Ngg[10];
		coincList[source_number].clear();
		
		for(int g0 = gammaList[source_number].size()-1; g0 >= 0 ; g0--){
				auto First = gammaList[source_number].at(g0);
				if( First.finalLvl == 0 ) continue;
				Ngg[0] = First.lvlPop * First.gammaBr;
				for(int g1 = g0-1; g1 >= 0; g1--){
						auto Second = gammaList[source_number].at(g1);
						if( First.finalLvl == Second.lvlEn ){
								Ngg[1] = Ngg[0] * Second.gammaBr;
								AddToCoincidenceList(source_number,{First},Second, Ngg[1] );
								if( Second.finalLvl == 0 ) continue;
								for(int g2 = g1-1; g2 >= 0; g2--){
										auto Third = gammaList[source_number].at(g2);
										if( Second.finalLvl == Third.lvlEn  ){
												Ngg[2] = Ngg[1] * Third.gammaBr;
												AddToCoincidenceList(source_number,{First,Second},Third, Ngg[2] );
												if( Third.finalLvl == 0 ) continue;
												for(int g3 = g2-1; g3 >= 0;g3--){
														auto Fourth = gammaList[source_number].at(g3);
														if( Third.finalLvl == Fourth.lvlEn ){
																Ngg[3] = Ngg[2] * Fourth.gammaBr;
																AddToCoincidenceList(source_number,{First,Second,Third},Fourth, Ngg[3] );
																if( Fourth.finalLvl == 0 ) continue;
																for(int g4 = g3 - 1; g4 >= 0; g4--){
																		auto Fifth = gammaList[source_number].at(g4);
																		if( Fourth.finalLvl == Fifth.lvlEn ){
																				Ngg[4] = Ngg[3] * Fifth.gammaBr;
																				AddToCoincidenceList(source_number,{First,Second,Third,Fourth},Fifth, Ngg[4] );
																				if( Fifth.finalLvl == 0 ) continue;
																				for(int g5 = g4-1; g5 >= 0; g5--){
																						auto Sixth = gammaList[source_number].at(g5);
																						if( Fifth.finalLvl == Sixth.lvlEn ){
																								Ngg[5] = Ngg[4] * Sixth.gammaBr;
																								AddToCoincidenceList(source_number,{First,Second,Third,Fourth,Fifth},Sixth, Ngg[5] );
																						} //end of g5 mapping
																				}//end of g5 loop
																		}//end of g4 mapping
																}//end of g4 loop
														}//end of g3 mapping
												}//end of g3 loop
										}//end of g2 map
								} //end of  g2 loop
						} //end of g1 match
				} //end of g1 loop
		}//end for g0 loop
		
}*/




/*
//this function gets a list of levels from the 'AssignedTransition' vector from the ReadInputData.cc file
//this function requires that ReadDecayScheme() function has been utilized
void GetLevelList(int source_number){

		if( gammaList[source_number].size() == 0 ){
				cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
				return;
		}

		double prev_energy = -1;
		int counter = 0;
		
		for(auto MyGamma : gammaList[source_number] ){
				if( MyGamma.lvlEn != prev_energy){
						MyLevels level = {counter, MyGamma.lvlEn, 0};
						levelList[source_number].push_back( level );
						counter++;
						prev_energy = MyGamma.lvlEn;
				}
		}
}

//this function is used to correct the 'AssignedTransition' list
//this correction is required if your final level energies do not exactly match the level energies in the 'LevelEnergy' vector
void FixTransitionList(int source_number){

		if( gammaList[source_number].size() == 0 ){
				cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
				return;
		}

		double min_diff;
		double diff;
		double newLvlEnergy;
		for(auto &MyGamma : gammaList[source_number]){
				min_diff = 10000.;
				for( auto level : levelList[source_number]){
						diff = TMath::Abs( MyGamma.finalLvl - level.lvlEn);
						if( diff < min_diff ){
								newLvlEnergy = level.lvlEn;
								min_diff = diff;
						}
				}
				MyGamma.finalLvl = newLvlEnergy;
		}

}


//this function calculates the level population of each state and the gamma-ray branching ratio of each transition
//this is function is required to properly calculate the expected gamma-gamma coincidence intensity
void CalcLevelFeedingAndGammaBR(int source_number){
	
		double sum_temp;
		for( auto &level : levelList[source_number]){
				sum_temp = 0;
				for(auto MyGamma : gammaList[source_number]){
						if( level.lvlEn == MyGamma.lvlEn ) sum_temp += MyGamma.gammaInt;
				}
				level.lvlPop = sum_temp;
		}
		
		//calculating gamma-ray branching ratios
		for(auto &MyGamma : gammaList[source_number]){
				if( MyGamma.lvlEn == 0) continue; //ground state has no decays so no branching ratio
				for( auto level : levelList[source_number]){
						if( level.lvlEn == MyGamma.lvlEn ){
								MyGamma.lvlPop = level.lvlPop;
								MyGamma.gammaBr = MyGamma.gammaInt / level.lvlPop;
						}
				}
		}
		
		//calculate level population
		for(auto &MyGamma : gammaList[source_number]){
				if( MyGamma.lvlEn == 0) continue; //ground state has no decays so no branching ratio
				for(auto otherGamma : gammaList[source_number]){
						if( otherGamma.finalLvl == MyGamma.lvlEn )
								MyGamma.lvlPop -= otherGamma.gammaInt;
								if( MyGamma.lvlPop < 0 ) MyGamma.lvlPop = 0;
				}
		}
}
*/
