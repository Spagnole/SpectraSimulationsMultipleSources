#include "../SimulateCoincidences.cc"

TH1D *hMultiplicity[NSources];
TH2D *hMultiVsEnergy[NSources];
double ggInt[10];
double gamma_multiplicity[100] = {0};
double NDecays;

void GetMultiplicity(int source_number, double GroundStateBranchRatio){

	cout << AssignedTransition[source_number].size() << endl;
	
	GetLevelList(source_number);
	FixTransitionList(source_number);
	CalcLevelFeedingAndGammaBR(source_number);

	

	hMultiplicity[source_number] = new TH1D( Form("hMult_source%d",source_number),  
					Form("hMult_source%d",source_number),
					100,-0.5,99.5);
	hMultiVsEnergy[source_number] = new TH2D( Form("hMultVsEnergy_source%d",source_number),  
					Form("hMultVsEnergy_source%d",source_number),
					100,-0.5,99.5,1000, -5, 9995);
		for(int i = AssignedTransition[source_number].size()-1; i >= 0 ; i--){
		
		NDecays = get<4>(AssignedTransition[source_number].at(i))*get<5>(AssignedTransition[source_number].at(i));
		if( get<2>(AssignedTransition[source_number].at(i))  == 0 ){
			gamma_multiplicity[1] += NDecays;
			hMultiVsEnergy[source_number]->Fill( 1, get<0>(AssignedTransition[source_number].at(i)), NDecays );
			continue;
		}
		for(int j = i-1; j >= 0 ; j--){
			if( get<2>(AssignedTransition[source_number].at(i)) == get<0>(AssignedTransition[source_number].at(j))){
				ggInt[2] = NDecays*get<5>(AssignedTransition[source_number].at(j));							
				if( get<2>(AssignedTransition[source_number].at(j))  == 0 ){ 
					gamma_multiplicity[2] += ggInt[2];
					hMultiVsEnergy[source_number]->Fill( 2, get<0>(AssignedTransition[source_number].at(i)), ggInt[2] );
					continue;
				}
				for(int k = j-1; k >= 0 ; k--){
					if( get<2>(AssignedTransition[source_number].at(j)) == get<0>(AssignedTransition[source_number].at(k))){
						ggInt[3] = ggInt[2]*get<5>(AssignedTransition[source_number].at(k));
						if( get<2>(AssignedTransition[source_number].at(k))  == 0 ){
							gamma_multiplicity[3] += ggInt[3];
							hMultiVsEnergy[source_number]->Fill( 3, get<0>(AssignedTransition[source_number].at(i)), ggInt[3] );
							continue;						
						}
						for(int l = k-1; l >= 0 ; l--){
							if( get<2>(AssignedTransition[source_number].at(k)) == get<0>(AssignedTransition[source_number].at(l))){
								ggInt[4] = ggInt[3]*get<5>(AssignedTransition[source_number].at(l));
								if( get<2>(AssignedTransition[source_number].at(l))  == 0 ){
									gamma_multiplicity[4] += ggInt[4];
									hMultiVsEnergy[source_number]->Fill( 4, get<0>(AssignedTransition[source_number].at(i)), ggInt[4] );
									continue;
								}
								for(int m = l-1; m >= 0 ; m--){
									if( get<2>(AssignedTransition[source_number].at(l)) == get<0>(AssignedTransition[source_number].at(m))){
										ggInt[5] = ggInt[4]*get<5>(AssignedTransition[source_number].at(m));
										if( get<2>(AssignedTransition[source_number].at(m))  == 0 ){
											gamma_multiplicity[5] += ggInt[5];
											hMultiVsEnergy[source_number]->Fill( 5, get<0>(AssignedTransition[source_number].at(i)), ggInt[5] );
											continue;
										}
										for(int n = m-1; n >= 0 ; n--){
											if( get<2>(AssignedTransition[source_number].at(m)) == get<0>(AssignedTransition[source_number].at(n))){
												ggInt[6] = ggInt[5]*get<5>(AssignedTransition[source_number].at(n));
												if( get<2>(AssignedTransition[source_number].at(n))  == 0 ){
													gamma_multiplicity[6] += ggInt[6];
													hMultiVsEnergy[source_number]->Fill( 6, get<0>(AssignedTransition[source_number].at(i)), ggInt[6] );
													continue;
												}
												for(int p = n-1; p >= 0 ; p--){
													if( get<2>(AssignedTransition[source_number].at(n)) == get<0>(AssignedTransition[source_number].at(p))){
														ggInt[7] = ggInt[6]*get<5>(AssignedTransition[source_number].at(p));
														if( get<2>(AssignedTransition[source_number].at(p)) == 0){
															gamma_multiplicity[7] += ggInt[7];
															hMultiVsEnergy[source_number]->Fill( 7, get<0>(AssignedTransition[source_number].at(i)), ggInt[7] );
															continue;
														}
														for(int q = p-1; q >= 0 ; q--){
															if( get<2>(AssignedTransition[source_number].at(p)) == get<0>(AssignedTransition[source_number].at(q))){
																ggInt[8] = ggInt[7]*get<5>(AssignedTransition[source_number].at(q));
																if( get<2>(AssignedTransition[source_number].at(q)) == 0){
																	gamma_multiplicity[8] += ggInt[8];
																	hMultiVsEnergy[source_number]->Fill( 8, get<0>(AssignedTransition[source_number].at(i)), ggInt[8] );
																	continue;
																}
																for(int r = q-1; r >= 0 ; r--){
																	if( get<2>(AssignedTransition[source_number].at(q)) == get<0>(AssignedTransition[source_number].at(r))){
																		ggInt[9] = ggInt[8]*get<5>(AssignedTransition[source_number].at(r));
																		if( get<2>(AssignedTransition[source_number].at(r)) == 0){
																			gamma_multiplicity[9] += ggInt[9];
																			hMultiVsEnergy[source_number]->Fill( 9, get<0>(AssignedTransition[source_number].at(i)), ggInt[9] );
																			continue;
																		}
																	}
																} //end of r loop
															}
														}	//end of q loop
													}
												} //end of p loop
											}
										} //end of n loop
									}
								} //end of m loop
							}
						} //end of l loop
					}
				} //end of k loop
			}
		} //end of j loop
	} //end of i loop
	
	//playing with matrix
	hMultiVsEnergy[source_number]->SetBinContent(1,1,hMultiVsEnergy[source_number]->Integral() / (1-GroundStateBranchRatio) );

	double total_int = 0;
	double NoGammasMult;
	double GammasMult = 1 - GroundStateBranchRatio;	
	for(int i = 1; i < 10; i++){
		total_int += gamma_multiplicity[i];
	}
	

	for(int i = 1; i < 10; i++){		
		hMultiplicity[source_number]->Fill( i,  GammasMult * gamma_multiplicity[i] / total_int  );
	}
	GammasMult = hMultiplicity[source_number]->Integral();
	hMultiplicity[source_number]->SetBinContent(1, 1-hMultiplicity[source_number]->Integral() );
	hMultiplicity[source_number]->Draw("hist");
	for(int i = 0; i < 10; i++){
		cout << "Normalised multiplicity " << i << " = "  <<  hMultiplicity[source_number]->GetBinContent(i+1) << endl;
	}
	cout << "Integral = " << hMultiplicity[source_number]->Integral() << endl;
}

vector<int> colors;
void make_colors_list( vector<int> color_num = {1,600,616,632,860+3,432,880,416+2,800,900,600-7,820-6,840+9,616-7,860+7,632-7,820,840} ){

	for(int i = 0; i < color_num.size(); i++){
		colors.push_back( color_num.at(i) );
	}
}

TH1D *hEnergyMult[10];
void DrawMultiplicities(int source_number = 0){

	make_colors_list();
	THStack *hs = new THStack("hs","Gamma-ray multiplicity");
	for(int i = 9; i >= 0; i--){
		hEnergyMult[i] = hMultiVsEnergy[source_number]->ProjectionY( Form("hMult%d_energy",i), i+1, i+1 );
		hEnergyMult[i]->SetTitle(Form("N_{#gamma} = %d",i));
		hEnergyMult[i]->SetLineColor( colors.at(i) );
		hEnergyMult[i]->SetFillColor( colors.at(i) );
		hEnergyMult[i]->SetFillStyle(3000);
		if( i == 0 ||  hEnergyMult[i]->Integral() == 0) continue;
		hs->Add(hEnergyMult[i]);
	}
	hs->Draw("hist");
	gPad->BuildLegend();
}
