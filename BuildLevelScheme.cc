#ifndef BUILD_LEVEL_INFO
#define BUILD_LEVEL_INFO

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

#endif
