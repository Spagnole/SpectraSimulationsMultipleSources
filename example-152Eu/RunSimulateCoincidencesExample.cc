#include "../SimulateCoincidences.cc"




void RunSim_Gate1087(){

	ReadDecayScheme("152Sm.dat","152Sm");
	ReadDecayScheme("152Gd.dat","152Gd");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate1087");	
	GetSpectrumBackground(hRealSpectra,50);
	
	for(int i = 0; i < used_sources; i++){
		GetLevelList(i);
		FixTransitionList(i);
		CalcLevelFeedingAndGammaBR(i);
		FindCoincidences(i);
		ReduceGammaGammaList(i);
		FillCoincMatrix(i);		
		GateOnSimMat(i, 1087-5,1087+4);
		FillEscPeakSpec(i);
	}
	BuildSimuledSpectra();
	
}

void RunSim_Gate244(){

	ReadDecayScheme("152Sm.dat","152Sm");
	ReadDecayScheme("152Gd.dat","152Gd");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate244");	
	GetSpectrumBackground(hRealSpectra,50);
	
	
	for(int i = 0; i < used_sources; i++){
		GetLevelList(i);
		FixTransitionList(i);
		CalcLevelFeedingAndGammaBR(i);
		FindCoincidences(i);
		ReduceGammaGammaList(i);
		FillCoincMatrix(i);		
		GateOnSimMat(i, 245-2,245+1);
		FillEscPeakSpec(i);
	}
	BuildSimuledSpectra();
	
}

void RunSim_Gate489(){

	ReadDecayScheme("152Sm.dat","152Sm");
	ReadDecayScheme("152Gd.dat","152Gd");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate489");	
	GetSpectrumBackground(hRealSpectra,50);
	for(int i = 0; i < used_sources; i++){
		GetLevelList(i);
		FixTransitionList(i);
		CalcLevelFeedingAndGammaBR(i);
		FindCoincidences(i);
		ReduceGammaGammaList(i);
		FillCoincMatrix(i);		
		GateOnSimMat(i, 489-3,489+2);
		FillEscPeakSpec(i);
	}
	BuildSimuledSpectra();

	
}

void RunSim_Gate676(){

	ReadDecayScheme("152Sm.dat","152Sm");
	ReadDecayScheme("152Gd.dat","152Gd");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate676");	
	GetSpectrumBackground(hRealSpectra,50);
	for(int i = 0; i < used_sources; i++){
		GetLevelList(i);
		FixTransitionList(i);
		CalcLevelFeedingAndGammaBR(i);
		FindCoincidences(i);
		ReduceGammaGammaList(i);
		FillCoincMatrix(i);		
		GateOnSimMat(i, 676-5,676+4);
		FillEscPeakSpec(i);
	}
	BuildSimuledSpectra();
	
}

void RunSim_Gate689(){

	ReadDecayScheme("152Sm.dat");
	ReadDecayScheme("152Gd.dat");
	GetEfficiency("MyExpEffnew.dat");
	GetPeakWidth("PeakWidths.dat");
	ReadEscapePeaks("EscapePeaks.dat");
	GetRealSpectra("ExampleFile.root","hggE_Gate689");	
	GetSpectrumBackground(hRealSpectra,50);
	
	for(int i = 0; i < used_sources; i++){
		GetLevelList(i);
		FixTransitionList(i);
		CalcLevelFeedingAndGammaBR(i);
		FindCoincidences(i);
		ReduceGammaGammaList(i);
		FillCoincMatrix(i);		
		GateOnSimMat(i, 689-3,689+2);
		FillEscPeakSpec(i);
	}
	BuildSimuledSpectra();
	
}

