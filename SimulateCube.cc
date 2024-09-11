#include "SimulateCoincidences.cc"
#include "ReadInputData.hh"
#include "WriteData.cc"
#include "BuildLevelScheme.cc"


TTree *gggCube[NSources];
TH2D *gated_matrix;
TH1D *gated_spectrum;

void BuildCube(int source_number = 0){
		double one,two,three;
		gggCube[source_number] = new TTree();
		gggCube[source_number]->SetName( Form("ggg_%d",source_number) );
		gggCube[source_number]->Branch("one",&one,"one/D");
		gggCube[source_number]->Branch("two",&two,"two/D");
		gggCube[source_number]->Branch("three",&three,"three/D");
		
		double counts;
		double sigma[3];
		for(auto ggg : cubeList[source_number]){
				counts = ( (NDet-2) /NDet)* ggg.triples * gEff->Eval(ggg.gammaOne)  * gEff->Eval(ggg.gammaTwo)  * gEff->Eval(ggg.gammaThree);
				sigma[0] = fWidth->Eval(ggg.gammaOne);
				sigma[1] = fWidth->Eval(ggg.gammaTwo);
				sigma[2] = fWidth->Eval(ggg.gammaThree);
				for(int i = 0; i < counts; i++){
						one = gRandom->Gaus(ggg.gammaOne,sigma[0]);
						two = gRandom->Gaus(ggg.gammaTwo,sigma[1]);
						three = gRandom->Gaus(ggg.gammaThree,sigma[2]);
						gggCube[source_number]->Fill();
				}
		}
		//gggCube[source_number]->SetDirectory(0);
}


void AddToCubeList(int source_number, vector<MyTransition> feederList, MyTransition Decay, double intensity) {
		//feederList.push_back(Decay);
		MyTransition gamma[2];
		Cube ggg;
		
		for (int i = 0; i < feederList.size(); i++) {
				gamma[0] = feederList.at(i);
				for (int j = i + 1; j < feederList.size(); j++) {
						gamma[1] = feederList.at(j);
						ggg = {gamma[0].Index, gamma[1].Index, Decay.Index, gamma[0].gammaEn, gamma[1].gammaEn, Decay.gammaEn, intensity, feederList};
						cubeList[source_number].push_back(ggg);
				}
		}
}


void FindTripleCoincidencesRec(int source_number, vector<MyTransition> &feederList, MyTransition currentTransition, double currentIntensity, int index) {
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
						AddToCubeList(source_number, feederList, nextTransition, nextIntensity);

						FindTripleCoincidencesRec(source_number, feederList, nextTransition, nextIntensity, i);
						feederList.pop_back(); // Backtrack
				}
		}
}


void FindTripleCoincidences(int source_number) {
		if (gammaList[source_number].size() == 0) {
				cout << "No data in decay list!\nHave you read in the list of gamma-rays?" << endl;
				return;
		}

		double Ngg[10];
		cubeList[source_number].clear();
		
		// Iterate through each transition as the starting point
		for (int g0 = gammaList[source_number].size() - 1; g0 >= 0; --g0) {
				MyTransition First = gammaList[source_number].at(g0);
				if (First.finalLvl == 0) continue;  // Skip transitions ending at level 0
				Ngg[0] = First.lvlPop * First.gammaBr;
				for(int g1 = g0-1; g1 >= 0; g1--){
						auto Second = gammaList[source_number].at(g1);
						if( First.finalLvl == Second.lvlEn ){
								Ngg[1] = Ngg[0] * Second.gammaBr;
								if( Second.finalLvl == 0 ) continue;
								
								// Initialize an empty feeder list
								vector<MyTransition> feederList;
								feederList.push_back(First);
								
								// Start recursive search with current transition as the initial feeder
								FindTripleCoincidencesRec(source_number, feederList, Second, Ngg[1], g1);
						}
				}
		}
}

void GateOnTree(TTree *tt, double first_low, double first_high, double second_low = 0, double second_high = 0 ){

		if( first_low >= first_high ){
				cout << "Bad gate!!! first_low > first_high??? Objects will be empty" <<endl;
				return;
		}
		string mat_name = Form("ggmat_gate%d",(int)(first_low + first_high)/2  );
		string mat_title = Form("gamma-gamma matrix gated a %d < E < %d", (int)first_low, (int)first_high );
		gated_matrix = new TH2D( mat_name.c_str() ,mat_title.c_str(), 1500, -0.5, 1499.5, 1500, -0.5, 1499.5 );
		
		string hist_name = Form("hSim_Gate_%d_Gate_%d",(int)(first_low + first_high)/2,(int)(second_low + second_high)/2);
		string hist_title = Form("Double Gated gamma spectrum, Gate: %d < E < %d, Gate: %d < E < %d",(int)first_low, (int)first_high, (int)second_low, (int)second_high );
		gated_spectrum = new TH1D(hist_name.c_str(),hist_title.c_str(),2500,-0.5,2499.5);
		double one, two, three;
		tt->SetBranchAddress("one", &one);
		tt->SetBranchAddress("two", &two);
		tt->SetBranchAddress("three", &three);
				
		for(int i = 0; i < tt->GetEntries(); i++){
				tt->GetEntry(i);
				if( one > first_low && one < first_high ){
						gated_matrix->Fill( two, three);
						gated_matrix->Fill( three, two);
						if( two > second_low && two < second_high )			gated_spectrum->Fill(three);
						if( three > second_low && three < second_high )	gated_spectrum->Fill(two);
				}
				if( two > first_low && two < first_high ){ //repeat of first if() statement with order changed
						gated_matrix->Fill( one, three);
						gated_matrix->Fill( three, one);
						if( one > second_low && one < second_high )			gated_spectrum->Fill(three);
						if( three > second_low && three < second_high )	gated_spectrum->Fill(one);
				}
				if( three > first_low && three < first_high ){ //repeat of first if() statement with order changed
						gated_matrix->Fill( one, two);
						gated_matrix->Fill( two, one);
						if( two > second_low && two < second_high )			gated_spectrum->Fill(one);
						if( one > second_low && one < second_high )			gated_spectrum->Fill(two);
				}
				
		}
}

void GateOnTree(double first_low, double first_high, double second_low = 0, double second_high = 0, int source_number = 0){
		GateOnTree( gggCube[source_number], first_low, first_high, second_low, second_high);
}

TFile *fWrite;
void WriteTree(string filename, int source_number = 0){
		fWrite = TFile::Open(filename.c_str(), "RECREATE");
		gggCube[source_number]->Write();
}

TFile *fTree;
void OpenTreeFile(string filename){
		fTree = TFile::Open(filename.c_str());
}

TFile *fSpectra;
void WriteSpectra(string filename, vector<TH1D*> hist_list = {0}, vector<TH2D*> mat_hist = {0}, TTree *tt = NULL){
		fSpectra = TFile::Open(filename.c_str(), "RECREATE");
		
		for(int i = 0; i < hist_list.size(); i++)	hist_list.at(i)->Write();
		for(int i = 0; i < mat_hist.size(); i++)	mat_hist.at(i)->Write();
		if( tt != NULL ) tt->Write();
}








//This function has been replaced by one which uses a recursive function to avoid the nested loops issue
//this function calculates the expected number of gamma-gamma-gamma coincidences
//this function uses a large number of nested loops
//if a gamma-gamma coincidence is separated by about 5 intermediate gamma-rays it will not be added to the list of coincidences
//if you require such coincidences ---> Add more nested loops
void OldFindTripleCoincidences(int source_number = 0){
		
		if( gammaList[source_number].size() == 0 ){
			cout << "No data in decay list!\nHave you read in list of gamma-rays?" << endl;
			return;
		}

		double Ngg[10];
		cubeList[source_number].clear();
		
		for(int g0 = gammaList[source_number].size()-1; g0 >= 0 ; g0--){
				auto First = gammaList[source_number].at(g0);
				if( First.finalLvl == 0 ) continue;
				Ngg[0] = First.lvlPop * First.gammaBr;
				for(int g1 = g0-1; g1 >= 0; g1--){
						auto Second = gammaList[source_number].at(g1);
						if( First.finalLvl == Second.lvlEn ){
								Ngg[1] = Ngg[0] * Second.gammaBr;
								if( Second.finalLvl == 0 ) continue;
								for(int g2 = g1-1; g2 >= 0; g2--){
										auto Third = gammaList[source_number].at(g2);
										if( Second.finalLvl == Third.lvlEn  ){
												Ngg[2] = Ngg[1] * Third.gammaBr;
												AddToCubeList(source_number,{First,Second},Third, Ngg[2] );
												if( Third.finalLvl == 0 ){
														continue;
												}
												for(int g3 = g2-1; g3 >= 0;g3--){
														auto Fourth = gammaList[source_number].at(g3);
														if( Third.finalLvl == Fourth.lvlEn ){
																Ngg[3] = Ngg[2] * Fourth.gammaBr;
																AddToCubeList(source_number,{First,Second,Third},Fourth, Ngg[3] );
																if( Fourth.finalLvl == 0 ){
																		continue;
																}
																for(int g4 = g3 - 1; g4 >= 0; g4--){
																		auto Fifth = gammaList[source_number].at(g4);
																		if( Fourth.finalLvl == Fifth.lvlEn ){
																				Ngg[4] = Ngg[3] * Fifth.gammaBr;
																				AddToCubeList(source_number,{First,Second,Third,Fourth},Fifth, Ngg[4] );
																				if( Fifth.finalLvl == 0 ){
																						continue;
																				}
																				for(int g5 = g4-1; g5 >= 0; g5--){
																						auto Sixth = gammaList[source_number].at(g5);
																						if( Fifth.finalLvl == Sixth.lvlEn ){
																								Ngg[5] = Ngg[4] * Sixth.gammaBr;
																								AddToCubeList(source_number,{First,Second,Third,Fourth,Fifth},Sixth, Ngg[5] );
																								if( Sixth.finalLvl == 0 ){
																										continue;
																								}
																								for(int g6 = g5-1; g6 >= 0; g6--){
																										auto Seventh = gammaList[source_number].at(g6);
																										if( Sixth.finalLvl == Seventh.lvlEn ){
																												Ngg[6] = Ngg[5] * Seventh.gammaBr;
																												AddToCubeList(source_number,{First,Second,Third,Fourth,Fifth,Sixth},Seventh, Ngg[6] );
																												if( Seventh.finalLvl == 0){
																														continue;
																												}
																												for(int g7 = g6-1; g7>=0; g7--){
																														auto Eighth = gammaList[source_number].at(g7);
																														if( Seventh.finalLvl == Eighth.lvlEn ){
																																Ngg[7] = Ngg[6] * Eighth.gammaBr;
																																AddToCubeList(source_number,{First,Second,Third,Fourth,Fifth,Sixth,Seventh},Eighth, Ngg[7] );
																																if( Eighth.finalLvl == 0){
																																		continue;
																																}
																																for(int g8 = g7-1; g8 >=0; g8--){
																																		auto Ninth = gammaList[source_number].at(g8);
																																		if( Eighth.finalLvl == Ninth.lvlEn ){
																																				Ngg[8] = Ngg[7] * Ninth.gammaBr;
																																				AddToCubeList(source_number,{First,Second,Third,Fourth,Fifth,Sixth,Seventh,Eighth},Ninth, Ngg[8] );
																																				if( Ninth.finalLvl == 0){
																																						continue;
																																				}
																																				for(int g9 = g8-1; g9 >= 0; g9--){
																																						auto Tenth = gammaList[source_number].at(g9);
																																						if( Ninth.finalLvl == Tenth.lvlEn ){
																																								Ngg[9] = Ngg[8] * Tenth.gammaBr;
																																								AddToCubeList(source_number,{First,Second,Third,Fourth,Fifth,Sixth,Seventh,Eighth,Ninth},Tenth, Ngg[9] );
																																								if( Tenth.finalLvl == 0){
																																										continue;
																																								}
																																						} //end of g9 mapping
																																				} //end of g9loop
																																		} //end of g9 mapping
																																} //end of g8 loop
																														} //end of g7 mapping
																												} //end of g7 loop
																										} //end of g6 mapping
																								} //end of g6 loop
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
		
}

// New function to check if the combination already exists in the cube list
//this was written for a previous version of AddToCube() which would try to insert repeat entries.
//current version does not have this issue.
bool CheckIfExistsInCubeList(int source_number, MyTransition gamma[3], const vector<MyTransition>& feederList) {
		Cube check_ggg;
		
		for (int check = 0; check < cubeList[source_number].size(); check++) { // Looping through existing cube list
				check_ggg = cubeList[source_number].at(check); // Get element from list
				if (gamma[0].Index == check_ggg.IndexOne && gamma[1].Index == check_ggg.IndexTwo &&
						gamma[2].Index == check_ggg.IndexThree && feederList.size() >= check_ggg.MyPath.size()) { // Check for same gammas
						
						// Compare feederList with the path in the cube list
						for (int second_check = 0; second_check < check_ggg.MyPath.size(); second_check++) {
								auto CheckMatch = feederList.at(second_check);
								auto Exists = check_ggg.MyPath.at(second_check);
								if (CheckMatch.Index != Exists.Index) {
										return false; // Mismatch found, continue the search
								}
						}
						return false; // Exact match found, skip adding
				}
		}
		
		return true; // New combination, proceed to add
}

