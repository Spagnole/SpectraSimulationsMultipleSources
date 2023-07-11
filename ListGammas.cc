#include "../ReadInputData.cc"

vector<tuple<double, double, double, double, double>> GammaList[NSources];


void MakeGammaList(int source_num){

	for(int i = 0; i < AssignedTransition[source_num].size(); i++){
		GammaList[source_num].push_back( make_tuple(get<1>(AssignedTransition[source_num].at(i)), get<0>(AssignedTransition[source_num].at(i)), get<2>(AssignedTransition[source_num].at(i)), get<3>(AssignedTransition[source_num].at(i)), get<4>(AssignedTransition[source_num].at(i)) ) );	
	}

	sort(GammaList[source_num].begin(), GammaList[source_num].end());

}


void PrintFullGammaList(int source_num){

		for(int i = 0; i < GammaList[source_num].size(); i++){
			cout << get<0>(GammaList[source_num].at(i)) << "\t" <<  get<1>(GammaList[source_num].at(i)) << "\t" <<   get<2>(GammaList[source_num].at(i)) << "\t" <<   get<3>(GammaList[source_num].at(i)) << "\t" <<   get<4>(GammaList[source_num].at(i)) << endl;	
	}
}

void SearchGammaList(int source_num, double low, double high){

		for(int i = 0; i < GammaList[source_num].size(); i++){
			if( get<0>(GammaList[source_num].at(i)) > low && get<0>(GammaList[source_num].at(i)) < high ){
				cout << get<0>(GammaList[source_num].at(i)) << "\t" <<  get<1>(GammaList[source_num].at(i)) << "\t" <<   get<2>(GammaList[source_num].at(i)) << "\t" <<   get<3>(GammaList[source_num].at(i)) << "\t" <<   get<4>(GammaList[source_num].at(i)) << endl;
		}
	}

}
